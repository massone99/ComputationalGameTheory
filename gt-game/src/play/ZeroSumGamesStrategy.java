package play;

import scpsolver.problems.LinearProgram;

import gametree.GameNode;
import gametree.GameTree;
import play.exception.InvalidStrategyException;
import scpsolver.constraints.Constraint;
import scpsolver.constraints.LinearBiggerThanEqualsConstraint;
import scpsolver.constraints.LinearEqualsConstraint;
import scpsolver.constraints.LinearSmallerThanEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;

import java.io.*;
import java.util.*;

public class ZeroSumGamesStrategy extends Strategy {

    public GameData computeGameData(GameNode rootNode) {
        int n1 = rootNode.numberOfChildren();
        int n2 = rootNode.getChildren().next().numberOfChildren();

        String[] labelsP1 = new String[n1];
        String[] labelsP2 = new String[n2];
        int[][] U1 = new int[n1][n2];
        int[][] U2 = new int[n1][n2];
        GameNode childNode1, childNode2;
        int j, i = 0;

        Iterator<GameNode> childrenNodes1 = rootNode.getChildren();

        while (childrenNodes1.hasNext()) {
            childNode1 = childrenNodes1.next();
            labelsP1[i] = childNode1.getLabel();
            j = 0;
            Iterator<GameNode> childrenNodes2 = childNode1.getChildren();
            while (childrenNodes2.hasNext()) {
                childNode2 = childrenNodes2.next();
                if (i == 0)
                    labelsP2[j] = childNode2.getLabel();
                U1[i][j] = childNode2.getPayoffP1();
                U2[i][j] = childNode2.getPayoffP2();
                j++;
            }
            i++;
        }

        GameData gameData = new GameData();
        gameData.U1 = U1;
        gameData.U2 = U2;
        gameData.labelsP1 = labelsP1;
        gameData.labelsP2 = labelsP2;

        return gameData;
    }

    public void execute() throws InterruptedException {

        while (!this.isTreeKnown()) {
            System.err.println("Waiting for game tree to become available.");
            Thread.sleep(1000);
        }

        serializeGameTree("EP5");

        while (true) {

            PlayStrategy myStrategy = this.getStrategyRequest();
            if (myStrategy == null) // Game was terminated by an outside event
                break;

            boolean playComplete = false; // Strategy was not valid

            while (!playComplete) {

                // Convert the game tree to a normal form game
                NormalFormGame nfg = gameTreeToNormalFormGame(tree.getRootNode());

                // Create linear programming problems for P1 and P2
                LP lpP1 = this.getLinearProblemP1(tree);
                LP lpP2 = this.getLinearProblemP2(tree);

                // Solve these problems to get the mixed strategies for P1 and P2
                LinearProgram lp1 = setLP1(lpP1.minFunCoeff, lpP1.constantTerms, lpP1.coefficientsMatrix, lpP1.lowerBounds);
                double[] mixedStrategyP1 = solveLP(lp1);
                LinearProgram lp2 = setLP1(lpP2.minFunCoeff, lpP2.constantTerms, lpP2.coefficientsMatrix, lpP2.lowerBounds);
                double[] mixedStrategyP2 = solveLP(lp2);
                showLP(lp1);
                showSolution(lp1, mixedStrategyP1);
                showLP(lp2);
                showSolution(lp2, mixedStrategyP2);


                // Printing Separator
                System.out.println("*******************************************************");

                // Init variables
                GameNode rootNode = tree.getRootNode();

                int n1 = rootNode.numberOfChildren();
                int n2 = rootNode.getChildren().next().numberOfChildren();

                GameData gameData = computeGameData(rootNode);

                String[] labelsP1 = gameData.labelsP1;
                String[] labelsP2 = gameData.labelsP2;
                int[][] U1 = gameData.U1;
                int[][] U2 = gameData.U2;

                System.out.println("ok");

                // Find payoffs for each possible pure strategy profile
                System.out.println("Assuming normal form game, see comment in the code.");
                Iterator<GameNode> childrenNodes1 = rootNode.getChildren();

                System.out.println("okok");

                Iterator<Integer> iterator = tree.getValidationSet().iterator();
                // Get iterator of moves already in the strategy
                Iterator<String> keys = myStrategy.keyIterator();

                double[] strategyP1 = new double[n1];
                double[] strategyP2 = new double[n2];

                for (int i = 0; i < labelsP1.length; i++) {
                    strategyP1[i] = i;
                    strategyP2[i] = i;
                }

                showStrategy(1, strategyP1, labelsP1);
                showStrategy(2, strategyP2, labelsP2);


                int i;
                for (i = 0; i < labelsP1.length; i++) {
                    myStrategy.put(labelsP1[i], mixedStrategyP1[i]);
                    myStrategy.put(labelsP2[i], mixedStrategyP2[i]);
                }

                // Print each strategy alongside his probability
                for (i = 0; i < labelsP1.length; i++) {
                    System.out.println(labelsP1[i] + "  " + mixedStrategyP1[i]);
                }

                for (i = 0; i < labelsP2.length; i++) {
                    System.out.println(labelsP2[i] + "  " + mixedStrategyP2[i]);
                }

                try {
                    this.provideStrategy(myStrategy);
                    playComplete = true;
                } catch (InvalidStrategyException e) {
                    System.err.println("Invalid strategy: " + e.getMessage());
                    e.printStackTrace(System.err);
                }
            }

        }
    }

    public String showLabel(String label) {
        return label.substring(label.lastIndexOf(':') + 1);
    }

    public void showStrategy(int P, double[] strategy, String[] labels) {
        System.out.println("Strategy Player " + P + ":");
        for (int i = 0; i < labels.length; i++) System.out.println("   " + strategy[i] + ":" + showLabel(labels[i]));
    }

    private void serializeGameTree(String gameName) {
        // Save the tree to a file
        try {
            FileOutputStream fileOut = new FileOutputStream(gameName + ".ser");
            ObjectOutputStream out = new ObjectOutputStream(fileOut);
            out.writeObject(tree);
            out.close();
            fileOut.close();
            System.out.printf("Serialized tree is saved in tree.ser");
        } catch (IOException i) {
            i.printStackTrace();
        }
    }

    public void printGame() {
        if (this.tree != null) {
            for (GameNode node : this.tree.getAllNodes()) {
                System.out.println(node);
            }
        } else {
            System.err.println("Game tree is not available.");
        }
    }

    private LP getLinearProblem(GameTree tree, boolean isPlayerOne) {
        double[][] equationsMatrix = isPlayerOne ? getEquationsMatrixP2(tree) : getEquationsMatrixP1(tree);
        double[][] equationsMatrixWithAddedConstant = addMinusOneColumnToMatrix(equationsMatrix);
        double[] minimizeFunctionCoefficients = getMinimizeFunctionCoefficients(equationsMatrixWithAddedConstant);
        int numConstantTerms = equationsMatrix[1].length + 1;
        double[] constantTerms = getConstantTerms(numConstantTerms);
        double[] lowerBounds = getLowerBounds(equationsMatrixWithAddedConstant);

        double[] constraintEquation = getConstraintOnVariablesEquation(equationsMatrixWithAddedConstant);
        double[][] coefficientsMatrix = addConstraintToEquationsMatrix(equationsMatrixWithAddedConstant, constraintEquation);
        return new LP(minimizeFunctionCoefficients, constantTerms, lowerBounds, coefficientsMatrix);
    }

    public LP getLinearProblemP1(GameTree tree) {
        return getLinearProblem(tree, true);
    }

    public LP getLinearProblemP2(GameTree tree) {
        return getLinearProblem(tree, false);
    }

    private double[][] addConstraintToEquationsMatrix(double[][] equationsMatrix, double[] constraintEquation) {
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

    private double[] getConstraintOnVariablesEquation(double[][] equationsMatrix) {
        int numVariables = equationsMatrix[0].length;
        double[] constraint = new double[numVariables];
        Arrays.fill(constraint, 1);

        // Set the last element to 0
        constraint[numVariables - 1] = 0;

        return constraint;
    }

    public void showLP(LinearProgram lp) {
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

    private double[] getLowerBounds(double[][] equationsMatrixP2) {
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

        return lowerBounds;
    }

    public GameTree deserializeGameTree() {
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

    public NormalFormGame gameTreeToNormalFormGame(GameNode rootNode) {
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

    public double[][] getPayoffP1ForAllStrategiesP2(GameNode rootNode) {
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

    public double[][] getPayoffP2ForAllStrategiesP1(GameNode rootNode) {
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

    public double[][] transposeMatrix(double[][] utilityMatrix) {
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

    public LinearProgram setLP1(double[] minOrMaxFun, double[] constantTerms, double[][] coeffiecientsOfVariables, double[] lowerBound) {
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

    public double[] solveLP(LinearProgram lp) {
        LinearProgramSolver solver = SolverFactory.newDefault();
        double[] x = solver.solve(lp);
        if (x == null) {
            throw new RuntimeException("No solution found");
        } else {
            return x;
        }
    }

    public void showSolution(LinearProgram lp, double[] x) {
        if (x == null) System.out.println("*********** NO SOLUTION FOUND ***********");
        else {
            System.out.println("*********** SOLUTION ***********");
            for (int i = 0; i < x.length; i++) System.out.println("x[" + i + "] = " + x[i]);
            System.out.println("f(x) = " + lp.evaluate(x));
        }
    }

    public double[][] utilityMatrixToEquations(double[][] utilityMatrix) {
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

    private double[] getConstantTerms(int numConstantTerms) {
        double[] constantTerms = new double[numConstantTerms];
        for (int i = 0; i < numConstantTerms - 1; i++) {
            constantTerms[i] = 0;
        }
        constantTerms[numConstantTerms - 1] = 1;
        return constantTerms;
    }

    private double[][] addMinusOneColumnToMatrix(double[][] equationsMatrixP2) {
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

    public double[] getMinimizeFunctionCoefficients(double[][] equationMatrix) {
        int length = equationMatrix[0].length;
        double[] result = new double[length];
        for (int i = 0; i < length - 1; i++) {
            result[i] = 0;
        }
        result[length - 1] = 1;
        return result;
    }

    private void printEquationsMatrix(double[][] equationsMatrix) {
        System.out.println("Equations matrix:");
        for (double[] row : equationsMatrix) {
            for (double value : row) {
                System.out.print(value + " ");
            }
            System.out.println();
        }
    }

    private double[][] getEquationsMatrixP2(GameTree tree) {
        double[][] payoffP2ForAllStrategiesP1 = getPayoffP2ForAllStrategiesP1(tree.getRootNode());
        return utilityMatrixToEquations(transposeMatrix(payoffP2ForAllStrategiesP1));
    }

    private double[][] getEquationsMatrixP1(GameTree tree) {
        double[][] payoffP1ForAllStrategiesP2 = getPayoffP1ForAllStrategiesP2(tree.getRootNode());
        return utilityMatrixToEquations(payoffP1ForAllStrategiesP2);
    }

    public boolean doProbabilitiesAddUpToOne(double[] probabilities) {
        double sum = 0.0;
        double epsilon = 0.00001;

        for (double probability : probabilities) {
            sum += probability;
        }

        return Math.abs(sum - 1.0) < epsilon;
    }

    public static class GameData {
        public int[][] U1;
        public int[][] U2;
        public String[] labelsP1;
        public String[] labelsP2;
    }

    public static class LP {
        public final double[] minFunCoeff;
        public final double[] constantTerms;
        public final double[] lowerBounds;
        public final double[][] coefficientsMatrix;

        public LP(double[] minFunCoeff, double[] constantTerms, double[] lowerBounds, double[][] coefficientsMatrix) {
            this.minFunCoeff = minFunCoeff;
            this.constantTerms = constantTerms;
            this.lowerBounds = lowerBounds;
            this.coefficientsMatrix = coefficientsMatrix;
        }
    }
}


/*
*
* package play;

import gametree.GameNode;
import gametree.GameTree;
import play.exception.InvalidStrategyException;
import scpsolver.constraints.Constraint;
import scpsolver.constraints.LinearBiggerThanEqualsConstraint;
import scpsolver.constraints.LinearEqualsConstraint;
import scpsolver.constraints.LinearSmallerThanEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;

import java.io.*;
import java.util.*;

public class MagicStrategy extends Strategy {

    public GameData computeGameData(GameNode rootNode) {
        int n1 = rootNode.numberOfChildren();
        int n2 = rootNode.getChildren().next().numberOfChildren();

        String[] labelsP1 = new String[n1];
        String[] labelsP2 = new String[n2];
        int[][] U1 = new int[n1][n2];
        int[][] U2 = new int[n1][n2];
        GameNode childNode1, childNode2;
        int j, i = 0;

        Iterator<GameNode> childrenNodes1 = rootNode.getChildren();

        while (childrenNodes1.hasNext()) {
            childNode1 = childrenNodes1.next();
            labelsP1[i] = childNode1.getLabel();
            j = 0;
            Iterator<GameNode> childrenNodes2 = childNode1.getChildren();
            while (childrenNodes2.hasNext()) {
                childNode2 = childrenNodes2.next();
                if (i == 0)
                    labelsP2[j] = childNode2.getLabel();
                U1[i][j] = childNode2.getPayoffP1();
                U2[i][j] = childNode2.getPayoffP2();
                j++;
            }
            i++;
        }

        GameData gameData = new GameData();
        gameData.U1 = U1;
        gameData.U2 = U2;
        gameData.labelsP1 = labelsP1;
        gameData.labelsP2 = labelsP2;

        return gameData;
    }

    public void execute() throws InterruptedException {

        while (!this.isTreeKnown()) {
            System.err.println("Waiting for game tree to become available.");
            Thread.sleep(1000);
        }

        //serializeGameTree();


        // This is where it should be solved i think?


        while (true) {

            PlayStrategy myStrategy = this.getStrategyRequest();
            if (myStrategy == null) // Game was terminated by an outside event
                break;

            boolean playComplete = false; // Strategy was not valid

            while (!playComplete) {

                // Convert the game tree to a normal form game
                NormalFormGame nfg = gameTreeToNormalFormGame(tree.getRootNode());

                // Create linear programming problems for P1 and P2
                LP lpP1 = this.getLinearProblemP1(tree);
                LP lpP2 = this.getLinearProblemP2(tree);

                // Solve these problems to get the mixed strategies for P1 and P2
                double[] mixedStrategyP1 = solveLP(setLP1(lpP1.minFunCoeff, lpP1.constantTerms, lpP1.coefficientsMatrix, lpP1.lowerBounds));
                double[] mixedStrategyP2 = solveLP(setLP1(lpP2.minFunCoeff, lpP2.constantTerms, lpP1.coefficientsMatrix, lpP2.lowerBounds));

                // Printing Separator
                System.out.println("*******************************************************");

                // Init variables
                GameNode rootNode = tree.getRootNode();

                int n1 = rootNode.numberOfChildren();
                int n2 = rootNode.getChildren().next().numberOfChildren();

                GameData gameData = computeGameData(rootNode);

                String[] labelsP1 = gameData.labelsP1;
                String[] labelsP2 = gameData.labelsP2;
                int[][] U1 = gameData.U1;
                int[][] U2 = gameData.U2;

                System.out.println("ok");

                // Find payoffs for each possible pure strategy profile
                System.out.println("Assuming normal form game, see comment in the code.");
                Iterator<GameNode> childrenNodes1 = rootNode.getChildren();

                System.out.println("okok");

                // Compute U1 and U2
//                while (childrenNodes1.hasNext()) {
//                    childNode1 = childrenNodes1.next();
//                    labelsP1[i] = childNode1.getLabel();
//                    j = 0;
//                    Iterator<GameNode> childrenNodes2 = childNode1.getChildren();
//                    while (childrenNodes2.hasNext()) {
//                        childNode2 = childrenNodes2.next();
//                        if (i == 0)
//                            labelsP2[j] = childNode2.getLabel(); // Assuming normal form game (perfect information not
//                        // supported)
//                        U1[i][j] = childNode2.getPayoffP1();
//                        U2[i][j] = childNode2.getPayoffP2();
//                        j++;
//                    }
//                    i++;
//                }

                // Get iterator of number of moves per turn

                Iterator<Integer> iterator = tree.getValidationSet().iterator();
                // Get iterator of moves already in the strategy
                Iterator<String> keys = myStrategy.keyIterator();

//                int i;
//
//                // Set random strategy
//                while (iterator.hasNext()) {
//                    Integer numStrategies = iterator.next();
//                    double[] moves = new double[numStrategies];
//                    double sum = 0;
//
//                    // For each move at this turn
//                    // Compute a random probability such that they sum up to 1
//                    for (i = 0; i < moves.length - 1; i++) {
//                        moves[i] = random.nextDouble();
//                        while (sum + moves[i] >= 1)
//                            moves[i] = random.nextDouble();
//                        sum = sum + moves[i];
//                    }
//                    moves[moves.length - 1] = ((double) 1) - sum;
//
//                    // Put in the strategy the probability for each action
//                    for (i = 0; i < moves.length; i++) {
//                        if (!keys.hasNext()) {
//                            System.err.println("PANIC: Strategy structure does not match the game.");
//                            return;
//                        }
//                        myStrategy.put(keys.next(), moves[i]);
//                    }
//                }

                double[] strategyP1 = new double[n1];
                double[] strategyP2 = new double[n2];

                showStrategy(1, strategyP1, labelsP1);
                showStrategy(2, strategyP2, labelsP2);


                int i;
                for (i = 0; i < labelsP1.length; i++) {
                    myStrategy.put(labelsP1[i], mixedStrategyP1[i]);
                    myStrategy.put(labelsP2[i], mixedStrategyP2[i]);
                }

                try {
                    this.provideStrategy(myStrategy);
                    playComplete = true;
                } catch (InvalidStrategyException e) {
                    System.err.println("Invalid strategy: " + e.getMessage());
                    e.printStackTrace(System.err);
                }
            }

        }
    }

    public String showLabel(String label) {
        return label.substring(label.lastIndexOf(':') + 1);
    }

    public void showStrategy(int P, double[] strategy, String[] labels) {
        System.out.println("Strategy Player " + P + ":");
        for (int i = 0; i < labels.length; i++) System.out.println("   " + strategy[i] + ":" + showLabel(labels[i]));
    }

    private void serializeGameTree() {
        // Save the tree to a file
        try {
            FileOutputStream fileOut = new FileOutputStream("tree.ser");
            ObjectOutputStream out = new ObjectOutputStream(fileOut);
            out.writeObject(tree);
            out.close();
            fileOut.close();
            System.out.printf("Serialized tree is saved in tree.ser");
        } catch (IOException i) {
            i.printStackTrace();
        }
    }

    public void printGame() {
        if (this.tree != null) {
            for (GameNode node : this.tree.getAllNodes()) {
                System.out.println(node);
            }
        } else {
            System.err.println("Game tree is not available.");
        }
    }

    private LP getLinearProblem(GameTree tree, boolean isPlayerOne) {
        double[][] equationsMatrix = isPlayerOne ? getEquationsMatrixP2(tree) : getEquationsMatrixP1(tree);
        double[][] equationsMatrixWithAddedConstant = addMinusOneColumnToMatrix(equationsMatrix);
        double[] minimizeFunctionCoefficients = getMinimizeFunctionCoefficients(equationsMatrixWithAddedConstant);
        int numConstantTerms = equationsMatrix[1].length + 1;
        double[] constantTerms = getConstantTerms(numConstantTerms);
        double[] lowerBounds = getLowerBounds(equationsMatrixWithAddedConstant);

        double[] constraintEquation = getConstraintOnVariablesEquation(equationsMatrixWithAddedConstant);
        double[][] coefficientsMatrix = addConstraintToEquationsMatrix(equationsMatrixWithAddedConstant, constraintEquation);
        return new LP(minimizeFunctionCoefficients, constantTerms, lowerBounds, coefficientsMatrix);
    }

    public LP getLinearProblemP1(GameTree tree) {
        return getLinearProblem(tree, true);
    }

    public LP getLinearProblemP2(GameTree tree) {
        return getLinearProblem(tree, false);
    }

    private double[][] addConstraintToEquationsMatrix(double[][] equationsMatrix, double[] constraintEquation) {
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

    private double[] getConstraintOnVariablesEquation(double[][] equationsMatrix) {
        int numVariables = equationsMatrix[0].length;
        double[] constraint = new double[numVariables];
        Arrays.fill(constraint, 1);

        // Set the last element to 0
        constraint[numVariables - 1] = 0;

        return constraint;
    }

    public void showLP(LinearProgram lp) {
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

    private double[] getLowerBounds(double[][] equationsMatrixP2) {
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

        return lowerBounds;
    }

    public GameTree deserializeGameTree() {
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

    public NormalFormGame gameTreeToNormalFormGame(GameNode rootNode) {
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

    public double[][] getPayoffP1ForAllStrategiesP2(GameNode rootNode) {
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

    public double[][] getPayoffP2ForAllStrategiesP1(GameNode rootNode) {
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

    public double[][] transposeMatrix(double[][] utilityMatrix) {
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

    public LinearProgram setLP1(double[] minOrMaxFun, double[] constantTerms, double[][] coeffiecientsOfVariables, double[] lowerBound) {
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

    public double[] solveLP(LinearProgram lp) {
        LinearProgramSolver solver = SolverFactory.newDefault();
        double[] x = solver.solve(lp);
        if (x == null) {
            throw new RuntimeException("No solution found");
        } else {
            return x;
        }
    }

    public void showSolution(LinearProgram lp, double[] x) {
        if (x == null) System.out.println("*********** NO SOLUTION FOUND ***********");
        else {
            System.out.println("*********** SOLUTION ***********");
            for (int i = 0; i < x.length; i++) System.out.println("x[" + i + "] = " + x[i]);
            System.out.println("f(x) = " + lp.evaluate(x));
        }
    }

    public double[][] utilityMatrixToEquations(double[][] utilityMatrix) {
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

    private double[] getConstantTerms(int numConstantTerms) {
        double[] constantTerms = new double[numConstantTerms];
        for (int i = 0; i < numConstantTerms - 1; i++) {
            constantTerms[i] = 0;
        }
        constantTerms[numConstantTerms - 1] = 1;
        return constantTerms;
    }

    private double[][] addMinusOneColumnToMatrix(double[][] equationsMatrixP2) {
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

    public double[] getMinimizeFunctionCoefficients(double[][] equationMatrix) {
        int length = equationMatrix[0].length;
        double[] result = new double[length];
        for (int i = 0; i < length - 1; i++) {
            result[i] = 0;
        }
        result[length - 1] = 1;
        return result;
    }

    private void printEquationsMatrix(double[][] equationsMatrix) {
        System.out.println("Equations matrix:");
        for (double[] row : equationsMatrix) {
            for (double value : row) {
                System.out.print(value + " ");
            }
            System.out.println();
        }
    }

    private double[][] getEquationsMatrixP2(GameTree tree) {
        double[][] payoffP2ForAllStrategiesP1 = getPayoffP2ForAllStrategiesP1(tree.getRootNode());
        return utilityMatrixToEquations(transposeMatrix(payoffP2ForAllStrategiesP1));
    }

    private double[][] getEquationsMatrixP1(GameTree tree) {
        double[][] payoffP1ForAllStrategiesP2 = getPayoffP1ForAllStrategiesP2(tree.getRootNode());
        return utilityMatrixToEquations(payoffP1ForAllStrategiesP2);
    }

    public boolean doProbabilitiesAddUpToOne(double[] probabilities) {
        double sum = 0.0;
        double epsilon = 0.00001;

        for (double probability : probabilities) {
            sum += probability;
        }

        return Math.abs(sum - 1.0) < epsilon;
    }

    public static class GameData {
        public int[][] U1;
        public int[][] U2;
        public String[] labelsP1;
        public String[] labelsP2;
    }

    public static class LP {
        public final double[] minFunCoeff;
        public final double[] constantTerms;
        public final double[] lowerBounds;
        public final double[][] coefficientsMatrix;

        public LP(double[] minFunCoeff, double[] constantTerms, double[] lowerBounds, double[][] coefficientsMatrix) {
            this.minFunCoeff = minFunCoeff;
            this.constantTerms = constantTerms;
            this.lowerBounds = lowerBounds;
            this.coefficientsMatrix = coefficientsMatrix;
        }
    }
}

* */