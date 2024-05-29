package play.tournament;

import gametree.GameNode;
import gametree.GameNodeDoesNotExistException;
import play.NormalFormGame;
import play.PlayStrategy;
import play.Strategy;
import play.exception.InvalidStrategyException;
import scpsolver.constraints.LinearSmallerThanEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;

import java.util.*;

public class T3_MixedFictitiousPlayStrategy extends Strategy {

    private Map<String, Integer> moveCountP2; // Count of player 2's moves
    private Map<String, Integer> moveCountP1; // Count of player 1's moves
    private int totalMovesP2; // Total moves count for player 2
    private int totalMovesP1; // Total moves count for player 1

    public static double[][] transposeMatrix(double[][] matrix) {
        int rows = matrix.length;
        int cols = matrix[0].length;
        double[][] transposed = new double[cols][rows];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                transposed[j][i] = matrix[i][j];
            }
        }
        return transposed;
    }

    @Override
    public void execute() throws InterruptedException {
        boolean isFirstRound;
        String[] actionsP2;
        String[] actionsP1;
        double[][] payoffMatrix2;
        while (!this.isTreeKnown()) {
            System.err.println("Waiting for game tree to become available!");
            Thread.sleep(1000);
        }

        GameNode rootNode = tree.getRootNode();
        GameInfo gameInfo = extractGameInfo(rootNode);
        // Payoff matrix for player 1
        double[][] payoffMatrix1 = gameInfo.payoffMatrix1;
        payoffMatrix2 = gameInfo.payoffMatrix2;
        actionsP1 = gameInfo.actionsP1;
        actionsP2 = gameInfo.actionsP2;

        moveCountP2 = new HashMap<>();
        moveCountP1 = new HashMap<>();
        for (String action : actionsP2) moveCountP2.put(action, 0);
        for (String action : actionsP1) moveCountP1.put(action, 0);
        totalMovesP2 = 0;
        totalMovesP1 = 0;
        isFirstRound = true;

        NormalFormGame normalFormGame = new NormalFormGame(payoffMatrix1, payoffMatrix2, actionsP1, actionsP2);
        normalFormGame.showGame();

        while (true) {
            PlayStrategy myStrategy = this.getStrategyRequest();
            if (myStrategy == null) break;

            boolean playComplete = false;
            while (!playComplete) {
                System.out.println("*******************************************************");

                GameNode lastMoveP1 = null;
                if (myStrategy.getFinalP1Node() != -1)
                    lastMoveP1 = this.tree.getNodeByIndex(myStrategy.getFinalP1Node());
                GameNode lastMoveP2 = null;
                if (myStrategy.getFinalP2Node() != -1)
                    lastMoveP2 = this.tree.getNodeByIndex(myStrategy.getFinalP2Node());

                double[] estimatedStrategyP2;
                double[] estimatedStrategyP1;

                if (lastMoveP1 != null && lastMoveP2 != null) {
                    try {
                        updateMoveCount(getReversePath(lastMoveP1), getReversePath(lastMoveP2));
                    } catch (GameNodeDoesNotExistException e) {
                        System.err.println("Game node does not exist: " + e.getMessage());
                        e.printStackTrace(System.err);
                    }

                    estimatedStrategyP2 = new double[actionsP2.length];
                    for (int i = 0; i < actionsP2.length; i++)
                        estimatedStrategyP2[i] = (double) moveCountP2.get(actionsP2[i]) / totalMovesP2;
                    estimatedStrategyP1 = new double[actionsP1.length];
                    for (int i = 0; i < actionsP1.length; i++)
                        estimatedStrategyP1[i] = (double) moveCountP1.get(actionsP1[i]) / totalMovesP1;

                    System.out.println("Estimated strategy for player 1: " + Arrays.toString(estimatedStrategyP1));
                    System.out.println("Estimated strategy for player 2: " + Arrays.toString(estimatedStrategyP2));

                } else {
                    estimatedStrategyP2 = new double[actionsP2.length];
                    Arrays.fill(estimatedStrategyP2, 1.0 / actionsP2.length);
                    estimatedStrategyP1 = new double[actionsP1.length];
                    Arrays.fill(estimatedStrategyP1, 1.0 / actionsP1.length);
                }

                double[] strategyP1;
                double[] strategyP2;

                if (isFirstRound) {
                    System.out.println("Using computeMixedStrategy for the first round.");
                    strategyP1 = computeMixedStrategy(payoffMatrix2);
                    strategyP2 = computeMixedStrategy(transposeMatrix(payoffMatrix1));
                    isFirstRound = false;
                } else {
                    System.out.println("Computing best response for player 1");
                    strategyP1 = computeBestResponse(payoffMatrix1, estimatedStrategyP2);
                    System.out.println("Computing best response for player 2");
                    strategyP2 = computeBestResponse(transposeMatrix(payoffMatrix2), estimatedStrategyP1);
                }

                showStrategy(1, strategyP1, actionsP1);
                showStrategy(2, strategyP2, actionsP2);

                for (int i = 0; i < actionsP1.length; i++)
                    myStrategy.put(actionsP1[i], strategyP1[i]);
                for (int i = 0; i < actionsP2.length; i++)
                    myStrategy.put(actionsP2[i], strategyP2[i]);

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

    private List<GameNode> getReversePath(GameNode node) {
        try {
            GameNode ancestor = node.getAncestor();
            List<GameNode> path = getReversePath(ancestor);
            path.add(node);
            return path;
        } catch (GameNodeDoesNotExistException e) {
            List<GameNode> path = new ArrayList<>();
            path.add(node);
            return path;
        }
    }

    private void updateMoveCount(List<GameNode> pathP1, List<GameNode> pathP2) throws GameNodeDoesNotExistException {
        List<String> opponentMoves = new ArrayList<>();

        for (GameNode node : pathP1) {
            if (node.isNature() || node.isRoot()) continue;
            if (node.getAncestor().isPlayer2()) {
                opponentMoves.add(node.getLabel());
            }
        }
        System.out.println("Last move for opponent of player 1: " + opponentMoves.get(opponentMoves.size() - 1));
        String lastMove = opponentMoves.get(opponentMoves.size() - 1);
        moveCountP2.put(lastMove, moveCountP2.get(lastMove) + 1);
        totalMovesP2++;

        for (GameNode node : pathP2) {
            if (node.isNature() || node.isRoot()) continue;
            if (node.getAncestor().isPlayer1()) {
                opponentMoves.add(node.getLabel());
            }
        }
        System.out.println("Last move for opponent of player 2: " + opponentMoves.get(opponentMoves.size() - 1));
        lastMove = opponentMoves.get(opponentMoves.size() - 1);
        moveCountP1.put(lastMove, moveCountP1.get(lastMove) + 1);
        totalMovesP1++;
    }

    private double[] computeBestResponse(double[][] payoffMatrix, double[] opponentStrategy) {
        double[] expectedPayoff = new double[payoffMatrix.length];
        for (int i = 0; i < payoffMatrix.length; i++) {
            for (int j = 0; j < payoffMatrix[i].length; j++) {
                expectedPayoff[i] += payoffMatrix[i][j] * opponentStrategy[j];
            }
        }
        System.out.println("Expected Payoff: " + Arrays.toString(expectedPayoff));

        double maxPayoff = Double.NEGATIVE_INFINITY;
        int bestAction = -1;
        for (int i = 0; i < expectedPayoff.length; i++) {
            if (expectedPayoff[i] > maxPayoff) {
                maxPayoff = expectedPayoff[i];
                bestAction = i;
            }
        }
        System.out.println("Best Action: " + bestAction);

        double[] bestResponse = new double[payoffMatrix.length];
        bestResponse[bestAction] = 1.0;
        return bestResponse;
    }

    private double[] computeMixedStrategy(double[][] utilities) {
        int n = utilities.length;
        int m = utilities[0].length;

        double[] coefficients = new double[n];
        Arrays.fill(coefficients, 1.0);

        double[][] constraints = new double[m + 1][n];
        double[] bounds = new double[m + 1];

        for (int j = 0; j < m; j++) {
            for (int i = 0; i < n; i++) {
                constraints[j][i] = utilities[i][j];
            }
            bounds[j] = 0;
        }

        for (int i = 0; i < n; i++) {
            constraints[m][i] = 1.0;
        }
        bounds[m] = 1.0;

        double[] lowerBounds = new double[n];
        Arrays.fill(lowerBounds, 0.0);

        LinearProgram lp = new LinearProgram(coefficients);
        lp.setMinProblem(false);
        for (int i = 0; i < bounds.length; i++) {
            lp.addConstraint(new LinearSmallerThanEqualsConstraint(constraints[i], bounds[i], "c" + i));
        }
        lp.setLowerbound(lowerBounds);

        double[] solution = solveLP(lp);
        if (solution == null) {
            throw new RuntimeException("No solution found for the linear program");
        }

        return solution;
    }

    public double[] solveLP(LinearProgram lp) {
        LinearProgramSolver solver = SolverFactory.newDefault();
        return solver.solve(lp);
    }

    public GameInfo extractGameInfo(GameNode rootNode) {
        int n1 = rootNode.numberOfChildren();
        int n2 = rootNode.getChildren().next().numberOfChildren();

        String[] actionsP1 = new String[n1];
        String[] actionsP2 = new String[n2];
        double[][] payoffMatrix1 = new double[n1][n2];
        double[][] payoffMatrix2 = new double[n1][n2];
        GameNode childNode1, childNode2;
        int j, i = 0;

        Iterator<GameNode> childrenNodes1 = rootNode.getChildren();

        while (childrenNodes1.hasNext()) {
            childNode1 = childrenNodes1.next();
            actionsP1[i] = childNode1.getLabel();
            j = 0;
            Iterator<GameNode> childrenNodes2 = childNode1.getChildren();
            while (childrenNodes2.hasNext()) {
                childNode2 = childrenNodes2.next();
                if (i == 0)
                    actionsP2[j] = childNode2.getLabel();
                payoffMatrix1[i][j] = childNode2.getPayoffP1();
                payoffMatrix2[i][j] = childNode2.getPayoffP2();
                j++;
            }
            i++;
        }

        GameInfo gameInfo = new GameInfo();
        gameInfo.payoffMatrix1 = payoffMatrix1;
        gameInfo.payoffMatrix2 = payoffMatrix2;
        gameInfo.actionsP1 = actionsP1;
        gameInfo.actionsP2 = actionsP2;

        return gameInfo;
    }

    public void showStrategy(int player, double[] strategy, String[] actions) {
        System.out.println("Strategy Player " + player + ":");
        for (int i = 0; i < actions.length; i++) {
            System.out.println("   " + strategy[i] + ": " + actions[i]);
        }
    }

    static class GameInfo {
        double[][] payoffMatrix1;
        double[][] payoffMatrix2;
        String[] actionsP1;
        String[] actionsP2;
    }
}
