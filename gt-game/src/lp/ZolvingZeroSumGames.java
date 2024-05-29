package lp;

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

public class ZolvingZeroSumGames {

    public static void main(String[] args) {
        GameTree tree = deserializeGameTree();
        assert tree != null;

        LinearProblem linearProblemP1 = getLinearProblem(tree, true);
        LinearProgram lp1 = setLP(linearProblemP1.minFunCoeff, linearProblemP1.constantTerms, linearProblemP1.coefficientsMatrix, linearProblemP1.lowerBounds);
        showLP(lp1);
        double[] x = solveLP(lp1);
        double[] mixedStrategy1 = Arrays.copyOfRange(x, 0, x.length - 1);
        System.out.println("Mixed strategy for player 1: " + Arrays.toString(mixedStrategy1));

        LinearProblem linearProblemP2 = getLinearProblem(tree, false);
        LinearProgram lp2 = setLP(linearProblemP2.minFunCoeff, linearProblemP2.constantTerms, linearProblemP2.coefficientsMatrix, linearProblemP2.lowerBounds);
        showLP(lp2);
        //till here everything right
        x = solveLP(lp2);
        double[] mixedStrategy2 = Arrays.copyOfRange(x, 0, x.length - 1);
        System.out.println("Mixed strategy for player 2: " + Arrays.toString(mixedStrategy2));
    }

    // must be right because equations are right
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

    // must be right because equations are right
    private static double[][] addConstraintToEquationsMatrix(double[][] equationsMatrix, double[] constraintEquation) {
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

    // must be right because equations are right
    private static double[] getConstraintOnVariablesEquation(double[][] equationsMatrix) {
        int numVariables = equationsMatrix[0].length;
        double[] constraint = new double[numVariables];
        Arrays.fill(constraint, 1);

        // Set the last element to 0
        constraint[numVariables - 1] = 0;

        return constraint;
    }

    // must be right because equations are right
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

    // must be right because equations are right
    private static double[] getLowerBounds(double[][] equationsMatrixP2) {
        int numConstantTerms = equationsMatrixP2[0].length;
        double[] lowerBounds = new double[numConstantTerms];

        double globalMin = Double.MAX_VALUE;
        for (double[] row : equationsMatrixP2) {
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

    // must be right because equations are right
    public static GameTree deserializeGameTree() {
        GameTree tree;
        try {
            FileInputStream fileIn = new FileInputStream("gt-game/tree.ser");
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

    // must be right because equations are right
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

    public static LinearProgram setLP(double[] minOrMaxFun, double[] constantTerms, double[][] coeffiecientsOfVariables, double[] lowerBound) {
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

    // must be right because equations are right
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

    // must be right because equations are right
    private static double[] getConstantTerms(int numConstantTerms) {
        double[] constantTerms = new double[numConstantTerms];
        for (int i = 0; i < numConstantTerms - 1; i++) {
            constantTerms[i] = 0;
        }
        constantTerms[numConstantTerms - 1] = 1;
        return constantTerms;
    }

    // must be right because equations are right
    private static double[][] addMinusOneColumnToMatrix(double[][] equationsMatrixP2) {
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

    // must be right because equations are right
    public static double[] getMinimizeFunctionCoefficients(double[][] equationMatrix) {
        int length = equationMatrix[0].length;
        double[] result = new double[length];
        for (int i = 0; i < length - 1; i++) {
            result[i] = 0;
        }
        result[length - 1] = 1;
        return result;
    }

    // must be right because equations are right
    private static double[][] getEquationsMatrixP2(GameTree tree) {
        double[][] payoffP2ForAllStrategiesP1 = getPayoffP2ForAllStrategiesP1(tree.getRootNode());
        return utilityMatrixToEquations(transposeMatrix(payoffP2ForAllStrategiesP1));
    }

    // must be right because equations are right
    private static double[][] getEquationsMatrixP1(GameTree tree) {
        double[][] payoffP1ForAllStrategiesP2 = getPayoffP1ForAllStrategiesP2(tree.getRootNode());
        return utilityMatrixToEquations(payoffP1ForAllStrategiesP2);
    }

    // must be right because equations are right

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

    // must be right because equations are right
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
