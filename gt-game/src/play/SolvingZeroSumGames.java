package play;

import gametree.GameNode;
import gametree.GameTree;
import scpsolver.constraints.Constraint;
import scpsolver.constraints.LinearBiggerThanEqualsConstraint;
import scpsolver.constraints.LinearEqualsConstraint;
import scpsolver.constraints.LinearSmallerThanEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.*;

public class SolvingZeroSumGames {

    public static void main(String[] args) {
        GameTree tree = deserializeGameTree("EP5");
        assert tree != null;
        NormalFormGame nfg = gameTreeToNormalFormGame(tree.getRootNode());
        nfg.showGame();

        LinearProblem linearProblemP1 = getLinearProblemP1(tree);

        LinearProgram lp1 = setLP1(
                linearProblemP1.minFunCoeff,
                linearProblemP1.constantTerms,
                linearProblemP1.coefficientsMatrix,
                linearProblemP1.lowerBounds);
        showLP(lp1);
        double[] x = solveLP(lp1);
        //showSolution(lp1, x);

        double[] mixedStrategy1 = Arrays.copyOfRange(x, 0, x.length - 1);
        System.out.println("Mixed strategy for player 1: " + Arrays.toString(mixedStrategy1));
        System.out.println();
        System.out.println("doProbabilitiesAddUpToOne(mixedStrategy1) = " + doProbabilitiesAddUpToOne(mixedStrategy1));

        LinearProblem linearProblemP2 = getLinearProblemP2(tree);

        LinearProgram lp2 = setLP1(linearProblemP2.minFunCoeff, linearProblemP2.constantTerms, linearProblemP2.coefficientsMatrix, linearProblemP2.lowerBounds);
        showLP(lp2);
        double[] x2 = solveLP(lp2);
        showSolution(lp2, x2);

        double[] mixedStrategy2 = Arrays.copyOfRange(x, 0, x.length - 1);
        System.out.println("Mixed strategy for player 2: " + Arrays.toString(mixedStrategy2));

        System.out.println();

        System.out.println("doProbabilitiesAddUpToOne(mixedStrategy1) = " + doProbabilitiesAddUpToOne(mixedStrategy2));
    }

    private static LinearProblem getLinearProblem(GameTree tree, boolean isPlayerOne) {
        double[][] equationsMatrix = isPlayerOne ? getEquationsMatrixP2(tree) : getEquationsMatrixP1(tree);
        double[][] equationsMatrixWithAddedConstant = addMinusOneColumnToMatrix(equationsMatrix);
        double[] minimizeFunctionCoefficients = getMinimizeFunctionCoefficients(equationsMatrixWithAddedConstant);
        int numConstantTerms = equationsMatrix[1].length + 1;
        double[] constantTerms = getConstantTerms(numConstantTerms);
        double[] lowerBounds = getLowerBounds(equationsMatrixWithAddedConstant);

        double[] constraintEquation = getConstraintOnVariablesEquation(equationsMatrixWithAddedConstant);
        double[][] coefficientsMatrix = addConstraintToEquationsMatrix(equationsMatrixWithAddedConstant, constraintEquation);
        return new LinearProblem(minimizeFunctionCoefficients, constantTerms, lowerBounds, coefficientsMatrix);
    }

    public static LinearProblem getLinearProblemP1(GameTree tree) {
        return getLinearProblem(tree, true);
    }

    public static LinearProblem getLinearProblemP2(GameTree tree) {
        return getLinearProblem(tree, false);
    }

    public static double[][] addConstraintToEquationsMatrix(double[][] equationsMatrix, double[] constraintEquation) {
        int numRows = equationsMatrix.length;
        int numCols = equationsMatrix[0].length;

        // Create a new matrix with an additional row
        double[][] newMatrix = new double[numRows + 1][numCols];

        // Copy the values from the old matrix to the new one
        for (int i = 0; i < numRows; i++) {
            System.arraycopy(equationsMatrix[i], 0, newMatrix[i], 0, numCols);
        }

        // Copy the constraint equation to the last row of the new matrix
        System.arraycopy(constraintEquation, 0, newMatrix[numRows], 0, numCols);

        return newMatrix;
    }

    public static double[] getConstraintOnVariablesEquation(double[][] equationsMatrix) {
        int numVariables = equationsMatrix[0].length;
        double[] constraint = new double[numVariables];
        Arrays.fill(constraint, 1);

        // Set the last element to 0
        constraint[numVariables - 1] = 0;

        return constraint;
    }

    public static void showLP(LinearProgram lp) {
        System.out.println("*********** LINEAR PROGRAMMING PROBLEM ***********");
        String fs;
        if (lp.isMinProblem()) System.out.print("  minimize: ");
        else System.out.print("  maximize: ");
        double[] cf = lp.getC();
        for (int i = 0; i < cf.length; i++)
            if (cf[i] != 0) {
                fs = String.format(Locale.US, "%+7.1f", cf[i]);
                System.out.print(fs + "*x[" + i + "]");
            }
        System.out.println();
        System.out.print("subject to: ");
        ArrayList<Constraint> lcstr = lp.getConstraints();
        double aij;
        double[] ci = null;
        String str = null;
        for (int i = 0; i < lcstr.size(); i++) {
            if (lcstr.get(i) instanceof LinearSmallerThanEqualsConstraint) {
                str = " <= ";
                ci = ((LinearSmallerThanEqualsConstraint) lcstr.get(i)).getC();
            }
            if (lcstr.get(i) instanceof LinearBiggerThanEqualsConstraint) {
                str = " >= ";
                ci = ((LinearBiggerThanEqualsConstraint) lcstr.get(i)).getC();
            }
            if (lcstr.get(i) instanceof LinearEqualsConstraint) {
                str = " == ";
                ci = ((LinearEqualsConstraint) lcstr.get(i)).getC();
            }
            str = str + String.format(Locale.US, "%6.1f", lcstr.get(i).getRHS());
            if (i != 0) System.out.print("            ");
            for (int j = 0; j < lp.getDimension(); j++) {
                assert ci != null;
                aij = ci[j];
                if (aij != 0) {
                    fs = String.format(Locale.US, "%+7.1f", aij);
                    System.out.print(fs + "*x[" + j + "]");
                } else System.out.print("            ");
            }
            System.out.println(str);
        }
    }

    public static double[] getLowerBounds(double[][] equationsMatrix) {
        int numConstantTerms = equationsMatrix[0].length;
        double[] lowerBounds = new double[numConstantTerms];

        double globalMin = Double.MAX_VALUE;
        for (double[] row : equationsMatrix) {
            for (double value : row) {
                if (value < globalMin) {
                    globalMin = value;
                }
            }
        }

        Arrays.fill(lowerBounds, 0);
        lowerBounds[numConstantTerms - 1] = globalMin;

        return lowerBounds;
    }

    public static GameTree deserializeGameTree(String gameName) {
        GameTree tree;
        try {
            FileInputStream fileIn = new FileInputStream("gt-game/" + gameName + ".ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            tree = (GameTree) in.readObject();
            in.close();
            fileIn.close();
        } catch (IOException i) {
            i.printStackTrace();
            return null;
        } catch (ClassNotFoundException c) {
            System.out.println("GameTree class not found");
            c.printStackTrace();
            return null;
        }
        return tree;
    }

    public static NormalFormGame gameTreeToNormalFormGame(GameNode rootNode) {
        int n1 = rootNode.numberOfChildren();
        int n2 = rootNode.getChildren().next().numberOfChildren();

        List<String> labelsP1 = new ArrayList<>();
        List<String> labelsP2 = new ArrayList<>();
        int[][] U1 = new int[n1][n2];
        int[][] U2 = new int[n1][n2];
        GameNode childNode1, childNode2;
        int j, i = 0;

        Iterator<GameNode> childrenNodes1 = rootNode.getChildren();

        while (childrenNodes1.hasNext()) {
            childNode1 = childrenNodes1.next();
            labelsP1.add(childNode1.getLabel());
            j = 0;
            Iterator<GameNode> childrenNodes2 = childNode1.getChildren();
            while (childrenNodes2.hasNext()) {
                childNode2 = childrenNodes2.next();
                if (i == 0) labelsP2.add(childNode2.getLabel());
                U1[i][j] = childNode2.getPayoffP1();
                U2[i][j] = childNode2.getPayoffP2();
                j++;
            }
            i++;
        }

        // Create the normal form game

        return new NormalFormGame(U1, U2, labelsP1.toArray(new String[0]), labelsP2.toArray(new String[0]));
    }

    public static double[][] getPayoffP1ForAllStrategiesP2(GameNode rootNode) {
        // Retrieve all the payoffs for player 1 for each strategy of player 2
        int n1 = rootNode.numberOfChildren();
        int n2 = rootNode.getChildren().next().numberOfChildren();

        // Create a 2D array to store the payoffs for player 1
        double[][] payoffsP1 = new double[n1][n2];

        // Get an iterator over the children of the root node, which represent the strategies of player 1
        Iterator<GameNode> strategiesP1 = rootNode.getChildren();
        int i = 0;
        while (strategiesP1.hasNext()) {
            GameNode strategyP1 = strategiesP1.next();

            // Iterate over the children of the current strategy of player 1
            Iterator<GameNode> childrenNodes = strategyP1.getChildren();
            int j = 0;
            while (childrenNodes.hasNext()) {
                GameNode childNode = childrenNodes.next();
                // Store the payoff for player 1 in the array
                payoffsP1[i][j] = childNode.getPayoffP1();
                j++;
            }
            i++;
        }

        return payoffsP1;
    }

    public static double[][] getPayoffP2ForAllStrategiesP1(GameNode rootNode) {
        // Retrieve all the payoffs for player 2 for each strategy of player 1
        int n1 = rootNode.numberOfChildren();
        int n2 = rootNode.getChildren().next().numberOfChildren();

        // Create a 2D array to store the payoffs for player 2
        double[][] payoffsP2 = new double[n1][n2];

        // Get an iterator over the children of the root node, which represent the strategies of player 1
        Iterator<GameNode> strategiesP1 = rootNode.getChildren();
        int i = 0;
        while (strategiesP1.hasNext()) {
            GameNode strategyP1 = strategiesP1.next();

            // Iterate over the children of the current strategy of player 1
            Iterator<GameNode> childrenNodes = strategyP1.getChildren();
            int j = 0;
            while (childrenNodes.hasNext()) {
                GameNode childNode = childrenNodes.next();
                // Store the payoff for player 2 in the array
                payoffsP2[i][j] = childNode.getPayoffP2();
                j++;
            }
            i++;
        }

        return payoffsP2;
    }

    public static double[][] transposeMatrix(double[][] utilityMatrix) {
        int numRows = utilityMatrix.length;
        int numCols = utilityMatrix[0].length;

        // Create a 2D array to store the transposed matrix
        double[][] transposedMatrix = new double[numCols][numRows];

        // Iterate over each row in the utility matrix
        for (int i = 0; i < numRows; i++) {
            // Iterate over each column in the utility matrix
            for (int j = 0; j < numCols; j++) {
                // Store the value from the utility matrix in the transposed matrix
                transposedMatrix[j][i] = utilityMatrix[i][j];
            }
        }

        return transposedMatrix;
    }

    public static LinearProgram setLP1(double[] minOrMaxFun, double[] constantTerms, double[][] coeffiecientsOfVariables, double[] lowerBound) {
        LinearProgram lp = new LinearProgram(minOrMaxFun);
        lp.setMinProblem(true);
        for (int i = 0; i < constantTerms.length; i++) {
            if (i == constantTerms.length - 1) {
                // Add the last constraint as an equality constraint
                lp.addConstraint(new LinearEqualsConstraint(coeffiecientsOfVariables[i], constantTerms[i], "coefficientsObjectiveFunction" + i));
            } else {
                // Add the other constraints as smaller than or equal to constraints
                lp.addConstraint(new LinearSmallerThanEqualsConstraint(coeffiecientsOfVariables[i], constantTerms[i], "coefficientsObjectiveFunction" + i));
            }
        }
        lp.setLowerbound(lowerBound);
        return lp;
    }

    public static double[] solveLP(LinearProgram lp) {
        LinearProgramSolver solver = SolverFactory.newDefault();
        double[] x = solver.solve(lp);
        if (x == null) {
            throw new RuntimeException("No solution found");
        } else {
            return x;
        }
    }

    public static void showSolution(LinearProgram lp, double[] x) {
        if (x == null) System.out.println("*********** NO SOLUTION FOUND ***********");
        else {
            System.out.println("*********** SOLUTION ***********");
            for (int i = 0; i < x.length; i++) System.out.println("x[" + i + "] = " + x[i]);
            System.out.println("f(x) = " + lp.evaluate(x));
        }
    }

    public static void showSolution(LinearProgram lp, double[] x, String[] l) {
        if (x == null)
            System.out.println("*********** NO SOLUTION FOUND ***********");
        else {
            System.out.println("*********** SOLUTION ***********");

            for (int i = 0; i < l.length; i++)
                // Print x[i] rounded at first 2 digits
                System.out.println(l[i] + " = " + String.format(Locale.US, "%.2f", x[i]));

            System.out.println("f(x) = " + lp.evaluate(x));

        }
    }

    public static double[][] utilityMatrixToEquations(double[][] utilityMatrix) {
        int numRows = utilityMatrix.length;
        int numCols = utilityMatrix[0].length;

        // Create a 2D array to store the matrix of equations for each column
        double[][] equationsMatrix = new double[numRows][numCols];

        // Iterate over each column in the utility matrix
        for (int j = 0; j < numCols; j++) {
            // Iterate over each row in the utility matrix
            for (int i = 0; i < numRows; i++) {
                // Store the value from the utility matrix in the equations matrix
                equationsMatrix[i][j] = utilityMatrix[i][j];
            }
        }

        return equationsMatrix;
    }

    public static double[] getConstantTerms(int numConstantTerms) {
        double[] constantTerms = new double[numConstantTerms];
        for (int i = 0; i < numConstantTerms - 1; i++) {
            constantTerms[i] = 0;
        }
        constantTerms[numConstantTerms - 1] = 1;
        return constantTerms;
    }

    public static double[][] addMinusOneColumnToMatrix(double[][] equationsMatrixP2) {
        int numRows = equationsMatrixP2.length;
        int numCols = equationsMatrixP2[0].length;

        // Create a new matrix with an additional column
        double[][] newMatrix = new double[numRows][numCols + 1];

        // Copy the values from the old matrix to the new one
        for (int i = 0; i < numRows; i++) {
            System.arraycopy(equationsMatrixP2[i], 0, newMatrix[i], 0, numCols);

            // Set the last element of each row to -1
            newMatrix[i][numCols] = -1;
        }

        return newMatrix;
    }

    public static double[] getMinimizeFunctionCoefficients(double[][] equationMatrix) {
        int length = equationMatrix[0].length;
        double[] result = new double[length];
        for (int i = 0; i < length - 1; i++) {
            result[i] = 0;
        }
        result[length - 1] = 1;
        return result;
    }

    private static void printEquationsMatrix(double[][] equationsMatrix) {
        System.out.println("Equations matrix:");
        for (double[] row : equationsMatrix) {
            for (double value : row) {
                System.out.print(value + " ");
            }
            System.out.println();
        }
    }

    public static double[][] getEquationsMatrixP2(GameTree tree) {
        double[][] payoffP2ForAllStrategiesP1 = getPayoffP2ForAllStrategiesP1(tree.getRootNode());
        return utilityMatrixToEquations(transposeMatrix(payoffP2ForAllStrategiesP1));
    }

    public static double[][] getEquationsMatrixP1(GameTree tree) {
        double[][] payoffP1ForAllStrategiesP2 = getPayoffP1ForAllStrategiesP2(tree.getRootNode());
        return utilityMatrixToEquations(payoffP1ForAllStrategiesP2);
    }

    public static boolean doProbabilitiesAddUpToOne(double[] probabilities) {
        double sum = 0.0;
        double epsilon = 0.00001;

        for (double probability : probabilities) {
            sum += probability;
        }

        return Math.abs(sum - 1.0) < epsilon;
    }

    public static class LinearProblem {
        public final double[] minFunCoeff;
        public final double[] constantTerms;
        public final double[] lowerBounds;
        public final double[][] coefficientsMatrix;

        public LinearProblem(double[] minFunCoeff, double[] constantTerms, double[] lowerBounds, double[][] coefficientsMatrix) {
            this.minFunCoeff = minFunCoeff;
            this.constantTerms = constantTerms;
            this.lowerBounds = lowerBounds;
            this.coefficientsMatrix = coefficientsMatrix;
        }
    }
}
