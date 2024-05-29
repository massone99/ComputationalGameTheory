package lp;

import gametree.GameTree;
import play.NormalFormGame;
import scpsolver.problems.LinearProgram;

import static play.SolvingZeroSumGames.*;

import java.util.Arrays;

public class MinMaxMaxMinApp {

    public static void main(String[] args) {
        // solveEM1typology();
        // solveEM4typology();
        solveEM7typology();
    }

    /**
     * Solve the EM7-9 typology of exercises
     * <p>
     * Main steps:
     * - creating a matrix of the utility values for player 2 with opposite sign
     * - solving the corresponding linear program
     */
    private static void solveEM7typology() {
        GameTree tree = deserializeGameTree("EM7");
        NormalFormGame nfg = gameTreeToNormalFormGame(tree.getRootNode());
        nfg.showGame();

        LinearProblem linearProblemP2 = getLinearProblemEM7(tree, true);
        LinearProgram lp2 = setLP1(
                linearProblemP2.minFunCoeff,
                linearProblemP2.constantTerms,
                linearProblemP2.coefficientsMatrix,
                linearProblemP2.lowerBounds);
        showLP(lp2);
        double[] x = solveLP(lp2);
        showSolution(lp2, x);

        double[] mixedStrategyP1 = Arrays.copyOfRange(x, 0, x.length - 1);
        System.out.println("Mixed Strategy P1: " + Arrays.toString(mixedStrategyP1));
    }

    /**
     * Solve the EM4-6 typology of exercises
     * <p>
     * Main steps:
     * - creating a matrix of the utility values for player 2 with opposite sign
     * - solving the corresponding linear program
     */
    private static void solveEM4typology() {
        GameTree tree = deserializeGameTree("EM4");
        NormalFormGame nfg = gameTreeToNormalFormGame(tree.getRootNode());
        nfg.showGame();

        LinearProblem linearProblemP2 = getInducedLinearProblemP1(tree);
        LinearProgram lp2 = setLP1(
                linearProblemP2.minFunCoeff,
                linearProblemP2.constantTerms,
                linearProblemP2.coefficientsMatrix,
                linearProblemP2.lowerBounds);
        showLP(lp2);
        double[] x = solveLP(lp2);
        showSolution(lp2, x);

        double[] mixedStrategyP1 = Arrays.copyOfRange(x, 0, x.length - 1);
        System.out.println("Mixed Strategy P1: " + Arrays.toString(mixedStrategyP1));
    }

    /**
     * Solve the EM1-3 typology of exercises
     * <p>
     * Main steps:
     * - creating a matrix of the utility values for player 1 with opposite sign
     * - transposing that matrix
     * - solving the corresponding linear program
     */
    private static void solveEM1typology() {
        GameTree tree = deserializeGameTree("EM3");
        NormalFormGame nfg = gameTreeToNormalFormGame(tree.getRootNode());
        nfg.showGame();

        // We need to create a new utility matrix by inverting the values of the payoff
        // for P1
        // Do not get confused by the name of the method, it is used to get the linear
        // problem for player 1
        LinearProblem linearProblemP1 = getInducedLinearProblemP2(tree);
        LinearProgram lp1 = setLP1(
                linearProblemP1.minFunCoeff,
                linearProblemP1.constantTerms,
                linearProblemP1.coefficientsMatrix,
                linearProblemP1.lowerBounds);
        showLP(lp1);
        double[] x = solveLP(lp1);
        showSolution(lp1, x);

        double[] mixedStrategyP1 = Arrays.copyOfRange(x, 0, x.length - 1);
        System.out.println("Mixed Strategy P1: " + Arrays.toString(mixedStrategyP1));

        double sum = 0;
        for (double value : mixedStrategyP1) {
            sum += value;
        }
        System.out.println("Sum of mixedStrategyP1: " + sum);
    }

    private static LinearProblem getInducedLinearProblemP2(GameTree tree) {
        return MinMaxMaxMinApp.getLinearProblem(tree, false);
    }

    private static LinearProblem getInducedLinearProblemP1(GameTree tree) {
        return MinMaxMaxMinApp.getLinearProblem(tree, true);
    }

    private static LinearProblem getLinearProblemEM7(GameTree tree, boolean isPlayerOne) {
        double[][] equationsMatrix = isPlayerOne ? getInducedUtilityMatrixP2EM7(tree) : getInducedUtilityMatrixP1(tree);
        // equationsMatrix = transposeMatrix(equationsMatrix);
        double[][] equationsMatrixWithAddedConstant = addMinusOneColumnToMatrix(equationsMatrix);
        double[] minimizeFunctionCoefficients = getMinimizeFunctionCoefficients(equationsMatrixWithAddedConstant);
        int numConstantTerms = equationsMatrix.length + 1;
        double[] constantTerms = getConstantTerms(numConstantTerms);
        double[] lowerBounds = getLowerBounds(equationsMatrixWithAddedConstant);

        double[] constraintEquation = getConstraintOnVariablesEquation(equationsMatrixWithAddedConstant);
        double[][] coefficientsMatrix = addConstraintToEquationsMatrix(equationsMatrixWithAddedConstant,
                constraintEquation);
        return new LinearProblem(minimizeFunctionCoefficients, constantTerms, lowerBounds, coefficientsMatrix);
    }

    private static LinearProblem getLinearProblem(GameTree tree, boolean isPlayerOne) {
        double[][] equationsMatrix = isPlayerOne ? getInducedUtilityMatrixP2(tree) : getInducedUtilityMatrixP1(tree);
        equationsMatrix = transposeMatrix(equationsMatrix);
        double[][] equationsMatrixWithAddedConstant = addMinusOneColumnToMatrix(equationsMatrix);
        double[] minimizeFunctionCoefficients = getMinimizeFunctionCoefficients(equationsMatrixWithAddedConstant);
        int numConstantTerms = equationsMatrix.length + 1;
        double[] constantTerms = getConstantTerms(numConstantTerms);
        double[] lowerBounds = getLowerBounds(equationsMatrixWithAddedConstant);

        double[] constraintEquation = getConstraintOnVariablesEquation(equationsMatrixWithAddedConstant);
        double[][] coefficientsMatrix = addConstraintToEquationsMatrix(equationsMatrixWithAddedConstant,
                constraintEquation);
        return new LinearProblem(minimizeFunctionCoefficients, constantTerms, lowerBounds, coefficientsMatrix);
    }

    private static double[][] getInducedUtilityMatrixP1(GameTree tree) {
        double[][] equationsMatrixP1 = getEquationsMatrixP1(tree);
        return negateMatrix(equationsMatrixP1);
    }

    private static double[][] getInducedUtilityMatrixP2(GameTree tree) {
        double[][] equationsMatrixP2 = getEquationsMatrixP2(tree);
        return negateMatrix(equationsMatrixP2);
    }

    private static double[][] getInducedUtilityMatrixP2EM7(GameTree tree) {
        double[][] equationsMatrixP2 = getEquationsMatrixP2(tree);
        return equationsMatrixP2;
    }

    private static double[][] negateMatrix(double[][] equationsMatrixP2) {
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
}
