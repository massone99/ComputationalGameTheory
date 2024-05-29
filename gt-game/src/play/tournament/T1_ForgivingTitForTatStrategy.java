package play.tournament;

import gametree.GameNode;
import gametree.GameNodeDoesNotExistException;
import play.PlayStrategy;
import play.Strategy;
import play.exception.InvalidStrategyException;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class T1_ForgivingTitForTatStrategy extends Strategy {

    private static final double FORGIVENESS_RATE = 0.2;
    private boolean defectStatusP1;
    private boolean defectStatusP2;
    private int roundCounter;
    private static final int MAX_ITERATIONS = 20; // Assume the game has 20 iterations

    private List<GameNode> traceBackPath(GameNode node) {
        try {
            GameNode ancestor = node.getAncestor();
            List<GameNode> path = traceBackPath(ancestor);
            path.add(node);
            return path;
        } catch (GameNodeDoesNotExistException e) {
            List<GameNode> path = new ArrayList<>();
            path.add(node);
            return path;
        }
    }

    private void evaluateStrategy(List<GameNode> pathP1, List<GameNode> pathP2) throws GameNodeDoesNotExistException {
        List<String> opponentMoves = new ArrayList<>();

        for (GameNode node : pathP1) {
            if (node.isNature() || node.isRoot()) continue;
            if (node.getAncestor().isPlayer2()) {
                opponentMoves.add(node.getLabel());
            }
        }

        if (opponentMoves.get(opponentMoves.size() - 1).contains("D")) {
            defectStatusP1 = Math.random() > FORGIVENESS_RATE;
        }

        opponentMoves.clear();
        for (GameNode node : pathP2) {
            if (node.isNature() || node.isRoot()) continue;
            if (node.getAncestor().isPlayer1()) {
                opponentMoves.add(node.getLabel());
            }
        }

        if (opponentMoves.get(opponentMoves.size() - 1).contains("D")) {
            defectStatusP2 = Math.random() > FORGIVENESS_RATE;
        }
    }

    private String extractLabel(String label) {
        return label.substring(label.lastIndexOf(':') + 1);
    }

    private void displayStrategy(int player, double[] strategy, String[] labels) {
        System.out.println("Strategy Player " + player + ":");
        for (int i = 0; i < labels.length; i++) {
            System.out.println("   " + strategy[i] + ":" + extractLabel(labels[i]));
        }
    }

    @Override
    public void execute() throws InterruptedException {

        defectStatusP1 = false;
        defectStatusP2 = false;
        roundCounter = 0;

        while (!this.isTreeKnown()) {
            System.err.println("Waiting for game tree to become available.");
            Thread.sleep(1000);
        }

        while (true) {

            roundCounter++;

            PlayStrategy currentStrategy = this.getStrategyRequest();
            if (currentStrategy == null) // Game was terminated by an outside event
                break;

            boolean roundComplete = false;

            while (!roundComplete) {

                System.out.println("*******************************************************");

                GameNode root = tree.getRootNode();
                int childrenCountP1 = root.numberOfChildren();
                int childrenCountP2 = root.getChildren().next().numberOfChildren();

                if (childrenCountP1 != 2 || childrenCountP2 != 2) {
                    System.err.println("ERROR: Game is not 2x2.");
                    return;
                }

                int i = 0;
                String[] movesP1 = new String[childrenCountP1];
                String[] movesP2 = new String[childrenCountP2];
                GameNode child1, child2;

                System.out.println("Assuming normal form game, see comment in the code.");
                Iterator<GameNode> iterator1 = root.getChildren();
                while (iterator1.hasNext()) {
                    child1 = iterator1.next();
                    movesP1[i] = child1.getLabel();
                    int j = 0;
                    Iterator<GameNode> iterator2 = child1.getChildren();
                    while (iterator2.hasNext()) {
                        child2 = iterator2.next();
                        if (i == 0)
                            movesP2[j] = child2.getLabel();
                        j++;
                    }
                    i++;
                }

                GameNode lastMoveP1 = null;
                GameNode lastMoveP2 = null;
                if (currentStrategy.getFinalP1Node() != -1) {
                    lastMoveP1 = this.tree.getNodeByIndex(currentStrategy.getFinalP1Node());
                }
                if (currentStrategy.getFinalP2Node() != -1) {
                    lastMoveP2 = this.tree.getNodeByIndex(currentStrategy.getFinalP2Node());
                }

                if (lastMoveP1 != null && lastMoveP2 != null) {
                    List<GameNode> pathP1 = traceBackPath(lastMoveP1);
                    List<GameNode> pathP2 = traceBackPath(lastMoveP2);
                    try {
                        evaluateStrategy(pathP1, pathP2);
                    } catch (GameNodeDoesNotExistException e) {
                        System.err.println("Game node does not exist: " + e.getMessage());
                        e.printStackTrace(System.err);
                    }
                }

                double[] strategyP1 = new double[childrenCountP1];
                double[] strategyP2 = new double[childrenCountP2];

                if (roundCounter == MAX_ITERATIONS) {
                    // Defect on the last iteration
                    strategyP1[0] = 0;
                    strategyP1[1] = 1;
                    strategyP2[0] = 0;
                    strategyP2[1] = 1;
                } else {
                    if (!defectStatusP1) {
                        strategyP1[0] = 1;
                        for (i = 1; i < movesP1.length; i++)
                            strategyP1[i] = 0;
                    } else {
                        for (i = 0; i < movesP1.length - 1; i++)
                            strategyP1[i] = 0;
                        strategyP1[movesP1.length - 1] = 1;
                    }

                    if (!defectStatusP2) {
                        strategyP2[0] = 1;
                        for (i = 1; i < movesP2.length; i++)
                            strategyP2[i] = 0;
                    } else {
                        for (i = 0; i < movesP2.length - 1; i++)
                            strategyP2[i] = 0;
                        strategyP2[movesP2.length - 1] = 1;
                    }
                }

                displayStrategy(1, strategyP1, movesP1);
                displayStrategy(2, strategyP2, movesP2);

                for (i = 0; i < movesP1.length; i++)
                    currentStrategy.put(movesP1[i], strategyP1[i]);
                for (i = 0; i < movesP2.length; i++)
                    currentStrategy.put(movesP2[i], strategyP2[i]);

                try {
                    this.provideStrategy(currentStrategy);
                    roundComplete = true;
                } catch (InvalidStrategyException e) {
                    System.err.println("Invalid strategy: " + e.getMessage());
                    e.printStackTrace(System.err);
                }
            }
        }
    }
}
