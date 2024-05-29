package play.tournament;

import gametree.GameNode;
import gametree.GameNodeDoesNotExistException;
import play.PlayStrategy;
import play.Strategy;
import play.exception.InvalidStrategyException;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class T2_AdaptiveCooperativeStrategy extends Strategy {

    private static final double TERMINATION_PROB = 0.05; // Likelihood that the current round is the final one
    private static final int MAX_ROUNDS = 100; // Upper limit of iterations
    private boolean lastMoveP1Defected = false; // Indicates if Player 1 defected in the previous round
    private boolean lastMoveP2Defected = false; // Indicates if Player 2 defected in the previous round
    private int roundCount = 0;

    @Override
    public void execute() throws InterruptedException {
        lastMoveP1Defected = false;
        lastMoveP2Defected = false;
        roundCount = 0;

        while (!this.isTreeKnown()) {
            System.err.println("Awaiting game tree availability.");
            Thread.sleep(1000);
        }

        while (true) {
            roundCount++;

            PlayStrategy currentStrategy = this.getStrategyRequest();
            if (currentStrategy == null) // Game terminated externally
                break;

            boolean isRoundComplete = false;

            while (!isRoundComplete) {

                System.out.println("*******************************************************");

                GameNode rootNode = tree.getRootNode();
                int childrenCountP1 = rootNode.numberOfChildren();
                int childrenCountP2 = rootNode.getChildren().next().numberOfChildren();

                if (childrenCountP1 != 2 || childrenCountP2 != 2) {
                    System.err.println("ERROR: Game is not a 2x2 game.");
                    return;
                }

                int i = 0;
                String[] actionsP1 = new String[childrenCountP1];
                String[] actionsP2 = new String[childrenCountP2];
                GameNode childNode1, childNode2;

                System.out.println("Assuming normal form game, see comment in the code.");
                Iterator<GameNode> childrenNodes1 = rootNode.getChildren();
                while (childrenNodes1.hasNext()) {
                    childNode1 = childrenNodes1.next();
                    actionsP1[i] = childNode1.getLabel();
                    int j = 0;
                    Iterator<GameNode> childrenNodes2 = childNode1.getChildren();
                    while (childrenNodes2.hasNext()) {
                        childNode2 = childrenNodes2.next();
                        if (i == 0)
                            actionsP2[j] = childNode2.getLabel();
                        j++;
                    }
                    i++;
                }

                GameNode finalMoveP1 = null;
                GameNode finalMoveP2 = null;
                if (currentStrategy.getFinalP1Node() != -1) {
                    finalMoveP1 = this.tree.getNodeByIndex(currentStrategy.getFinalP1Node());
                }
                if (currentStrategy.getFinalP2Node() != -1) {
                    finalMoveP2 = this.tree.getNodeByIndex(currentStrategy.getFinalP2Node());
                }

                if (finalMoveP1 != null && finalMoveP2 != null) {
                    List<GameNode> pathP1 = traceReversePath(finalMoveP1);
                    List<GameNode> pathP2 = traceReversePath(finalMoveP2);
                    try {
                        updateDefectionStatus(pathP1, pathP2);
                    } catch (GameNodeDoesNotExistException e) {
                        System.err.println("Game node not found: " + e.getMessage());
                        e.printStackTrace(System.err);
                    }
                }

                double[] probabilitiesP1 = new double[childrenCountP1];
                double[] probabilitiesP2 = new double[childrenCountP2];

                // Update strategy based on round number and continuation probability
                double continuationLikelihood = Math.pow(1 - TERMINATION_PROB, roundCount);
                if (roundCount >= MAX_ROUNDS || continuationLikelihood < 0.5) {
                    // Higher chance of game ending soon - favor defection
                    probabilitiesP1[0] = 0;
                    probabilitiesP1[1] = 1;
                    probabilitiesP2[0] = 0;
                    probabilitiesP2[1] = 1;
                } else {
                    // More rounds expected - favor cooperation
                    if (!lastMoveP1Defected) {
                        probabilitiesP1[0] = 1;
                        for (i = 1; i < actionsP1.length; i++)
                            probabilitiesP1[i] = 0;
                    } else {
                        for (i = 0; i < actionsP1.length - 1; i++)
                            probabilitiesP1[i] = 0;
                        probabilitiesP1[actionsP1.length - 1] = 1;
                    }

                    if (!lastMoveP2Defected) {
                        probabilitiesP2[0] = 1;
                        for (i = 1; i < actionsP2.length; i++)
                            probabilitiesP2[i] = 0;
                    } else {
                        for (i = 0; i < actionsP2.length - 1; i++)
                            probabilitiesP2[i] = 0;
                        probabilitiesP2[actionsP2.length - 1] = 1;
                    }
                }

                displayStrategy(1, probabilitiesP1, actionsP1);
                displayStrategy(2, probabilitiesP2, actionsP2);

                // Assuming it's a 2x2 game
                for (i = 0; i < actionsP1.length; i++)
                    currentStrategy.put(actionsP1[i], probabilitiesP1[i]);
                for (i = 0; i < actionsP2.length; i++)
                    currentStrategy.put(actionsP2[i], probabilitiesP2[i]);

                // Submit strategy
                try {
                    this.provideStrategy(currentStrategy);
                    isRoundComplete = true;
                } catch (InvalidStrategyException e) {
                    System.err.println("Invalid strategy: " + e.getMessage());
                    e.printStackTrace(System.err);
                }

            }
        }
    }

    private List<GameNode> traceReversePath(GameNode currentNode) {
        try {
            GameNode ancestor = currentNode.getAncestor();
            List<GameNode> path = traceReversePath(ancestor);
            path.add(currentNode);
            return path;
        } catch (GameNodeDoesNotExistException e) {
            List<GameNode> path = new ArrayList<>();
            path.add(currentNode);
            return path;
        }
    }

    private void updateDefectionStatus(List<GameNode> pathP1, List<GameNode> pathP2) throws GameNodeDoesNotExistException {
        // Check the last move of opponents and update defection status
        if (pathP1 != null && pathP1.size() > 1) {
            GameNode lastMoveP1 = pathP1.get(pathP1.size() - 2); // Last move made by P1
            lastMoveP1Defected = lastMoveP1.getLabel().equals("Defect");
        }

        if (pathP2 != null && pathP2.size() > 1) {
            GameNode lastMoveP2 = pathP2.get(pathP2.size() - 2); // Last move made by P2
            lastMoveP2Defected = lastMoveP2.getLabel().equals("Defect");
        }
    }

    private void displayStrategy(int player, double[] strategy, String[] labels) {
        System.out.println("Strategy for Player " + player + ":");
        for (int i = 0; i < strategy.length; i++) {
            System.out.println(labels[i] + ": " + strategy[i]);
        }
    }

}
