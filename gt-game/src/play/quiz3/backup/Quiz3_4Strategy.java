package play.quiz3.backup;


import gametree.GameNode;
import gametree.GameTree;
import play.GameData;
import play.NormalFormGame;
import play.PlayStrategy;
import play.Strategy;
import play.exception.InvalidStrategyException;
import scpsolver.problems.LinearProgram;

import static play.SolvingZeroSumGames.*;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.Iterator;

public class Quiz3_4Strategy extends Strategy {
    private static double[][] getInducedUtilityMatrixP2EM7(GameTree tree) {
        System.out.println("**** USING QUIZ 3 CODE FOR QUIZ4 **** getInducedUtilityMatrixP2EM7");
        double[][] equationsMatrixP2 = getEquationsMatrixP2(tree);
        return equationsMatrixP2;
    }

    public void serializeGameTree(GameTree tree, String filePath) {
        try {
            FileOutputStream fileOut = new FileOutputStream(filePath);
            ObjectOutputStream out = new ObjectOutputStream(fileOut);
            out.writeObject(tree);
            out.close();
            fileOut.close();
        } catch (IOException i) {
            i.printStackTrace();
        }
    }

    public void execute() throws InterruptedException {

        while (!this.isTreeKnown()) {
            System.err.println("Waiting for game tree to become available.");
            Thread.sleep(1000);
        }

        //serializeGameTree(this.tree, "EM10.ser");
        NormalFormGame nfg = gameTreeToNormalFormGame(tree.getRootNode());
        nfg.showGame();

        while (true) {

            PlayStrategy myStrategy = this.getStrategyRequest();
            if (myStrategy == null) // Game was terminated by an outside event
                break;

            boolean playComplete = false; // Strategy was not valid

            while (!playComplete) {

                //NormalFormGame nfg = gameTreeToNormalFormGame(tree.getRootNode());
                //nfg.showGame();

                LinearProblem linearProblemP2 = getLinearProblemEM10(tree, false);
                LinearProgram lp2 = setLP1(linearProblemP2.minFunCoeff, linearProblemP2.constantTerms, linearProblemP2.coefficientsMatrix, linearProblemP2.lowerBounds);
                showLP(lp2);
                double[] x = solveLP(lp2);
                showSolution(lp2, x);

                double[] mixedStrategyP1 = Arrays.copyOfRange(x, 0, x.length - 1);
                System.out.println("Mixed Strategy P1: " + Arrays.toString(mixedStrategyP1));

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

                // Get iterator of number of moves per turn
                Iterator<Integer> iterator = tree.getValidationSet().iterator();
                // Get iterator of strategies
                Iterator<String> keys = myStrategy.keyIterator();

                // Set strategy
                while (iterator.hasNext()) {
                    double[] moves = new double[iterator.next()];

                    // Set the probability of the first label to 1 and all the others to 0
                    for (int i = 0; i < moves.length; i++) {
                        moves[i] = (i == 0) ? 1.0 : 0.0;
                    }

                    // Put in the strategy the probability for each action
                    for (double move : moves) {
                        if (!keys.hasNext()) {
                            System.err.println("PANIC: Strategy structure does not match the game.");
                            break;
                        }
                        myStrategy.put(keys.next(), move);
                    }
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

    private LinearProblem getLinearProblemEM10(GameTree tree, boolean isPlayerOne) {
        double[][] equationsMatrix = isPlayerOne ? getInducedUtilityMatrixP2EM7(tree) : getInducedUtilityMatrixP1EM10(tree);
        double[][] equationsMatrixWithAddedConstant = addMinusOneColumnToMatrix(equationsMatrix);
        double[] minimizeFunctionCoefficients = getMinimizeFunctionCoefficients(equationsMatrixWithAddedConstant);
        int numConstantTerms = equationsMatrix.length + 1;
        double[] constantTerms = getConstantTerms(numConstantTerms);
        double[] lowerBounds = getLowerBounds(equationsMatrixWithAddedConstant);

        double[] constraintEquation = getConstraintOnVariablesEquation(equationsMatrixWithAddedConstant);
        double[][] coefficientsMatrix = addConstraintToEquationsMatrix(equationsMatrixWithAddedConstant, constraintEquation);
        return new LinearProblem(minimizeFunctionCoefficients, constantTerms, lowerBounds, coefficientsMatrix);
    }

    private double[][] getInducedUtilityMatrixP1EM10(GameTree tree) {
        double[][] equationsMatrixP1 = getEquationsMatrixP1(tree);
        return equationsMatrixP1;
    }

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
                if (i == 0) labelsP2[j] = childNode2.getLabel();
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

    public String showLabel(String label) {
        return label.substring(label.lastIndexOf(':') + 1);
    }

    public void showStrategy(int P, double[] strategy, String[] labels) {
        System.out.println("Strategy Player " + P + ":");
        for (int i = 0; i < labels.length; i++) System.out.println("   " + strategy[i] + ":" + showLabel(labels[i]));
    }

    public class LinearProblem {
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

