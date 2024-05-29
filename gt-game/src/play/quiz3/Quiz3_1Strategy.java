package play.quiz3;

import play.GameData;
import play.NormalFormGame;
import play.PlayStrategy;
import play.Strategy;
import scpsolver.problems.LinearProgram;
import gametree.GameNode;
import gametree.GameTree;
import play.exception.InvalidStrategyException;

import static play.SolvingZeroSumGames.*;

import java.util.*;

public class Quiz3_1Strategy extends Strategy {

    public void execute() throws InterruptedException {

        while (!this.isTreeKnown()) {
            System.err.println("Waiting for game tree to become available.");
            Thread.sleep(1000);
        }

        NormalFormGame nfg = gameTreeToNormalFormGame(tree.getRootNode());
        nfg.showGame();

        while (true) {

            PlayStrategy myStrategy = this.getStrategyRequest();
            if (myStrategy == null) // Game was terminated by an outside event
                break;

            boolean playComplete = false; // Strategy was not valid

            while (!playComplete) {

                // Init variables
                int i = 0;
                GameNode rootNode = tree.getRootNode();
                GameNode childNode1;
                int n1 = rootNode.numberOfChildren();
                String[] labels = new String[n1];
                Iterator<GameNode> childrenNodes1 = rootNode.getChildren();
                while (childrenNodes1.hasNext()) {
                    childNode1 = childrenNodes1.next();
                    labels[i] = childNode1.getLabel();
                    i++;
                }

                // We need to create a new utility matrix by inverting the values of the payoff for P1
                LinearProblem linearProblemP1 = getInducedLinearProblemP2(tree);
                LinearProgram lp1 = setLP1(linearProblemP1.minFunCoeff, linearProblemP1.constantTerms, linearProblemP1.coefficientsMatrix, linearProblemP1.lowerBounds);
                showLP(lp1);
                double[] x = solveLP(lp1);

                showSolution(lp1, x, labels);

                // Printing Separator
                System.out.println("*******************************************************");

                // Get iterator of number of moves per turn
                Iterator<Integer> iterator = tree.getValidationSet().iterator();
                // Get iterator of strategies
                Iterator<String> keys = myStrategy.keyIterator();

                // Set strategy
                while (iterator.hasNext()) {
                    double[] moves = new double[iterator.next()];

                    // Set the probability of the first label to 1 and all the others to 0
                    for (i = 0; i < moves.length; i++)
                        moves[i] = (i == 0) ? 1.0 : 0.0;

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


    private LinearProblem getInducedLinearProblemP2(GameTree tree) {
        return getLinearProblem(tree, false);
    }

    private LinearProblem getInducedLinearProblemP1(GameTree tree) {
        return getLinearProblem(tree, true);
    }

    private LinearProblem getLinearProblem(GameTree tree, boolean isPlayerOne) {
        double[][] equationsMatrix = isPlayerOne ? getInducedUtilityMatrixP2(tree) : getInducedUtilityMatrixP1(tree);
        equationsMatrix = transposeMatrix(equationsMatrix);
        double[][] equationsMatrixWithAddedConstant = addMinusOneColumnToMatrix(equationsMatrix);
        double[] minimizeFunctionCoefficients = getMinimizeFunctionCoefficients(equationsMatrixWithAddedConstant);
        int numConstantTerms = equationsMatrix.length + 1;
        double[] constantTerms = getConstantTerms(numConstantTerms);
        double[] lowerBounds = getLowerBounds(equationsMatrixWithAddedConstant);

        double[] constraintEquation = getConstraintOnVariablesEquation(equationsMatrixWithAddedConstant);
        double[][] coefficientsMatrix = addConstraintToEquationsMatrix(equationsMatrixWithAddedConstant, constraintEquation);
        return new LinearProblem(minimizeFunctionCoefficients, constantTerms, lowerBounds, coefficientsMatrix);
    }

    private double[][] getInducedUtilityMatrixP1(GameTree tree) {
        double[][] equationsMatrixP1 = getEquationsMatrixP1(tree);
        return negateMatrix(equationsMatrixP1);
    }

    private double[][] getInducedUtilityMatrixP2(GameTree tree) {
        double[][] equationsMatrixP2 = getEquationsMatrixP2(tree);
        return negateMatrix(equationsMatrixP2);
    }

    private double[][] negateMatrix(double[][] equationsMatrixP2) {
        double[][] result = new double[equationsMatrixP2.length][equationsMatrixP2[0].length];
        for (int i = 0; i < equationsMatrixP2.length; i++) {
            // Iterate over each column in the row
            for (int j = 0; j < equationsMatrixP2[i].length; j++) {
                // Multiply the element by -1 to negate it
                result[i][j] = equationsMatrixP2[i][j] * -1;
            }
        }
        return result;
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
