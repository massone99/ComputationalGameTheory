package play;

import gametree.GameNode;
import play.exception.InvalidStrategyException;
import scpsolver.constraints.LinearSmallerThanEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;

import java.util.Arrays;
import java.util.Iterator;

public class TournamentMixedNEStrategy extends Strategy {

    @Override
    public void execute() throws InterruptedException {
        // Wait until the game tree is available
        while (!this.isTreeKnown()) {
            System.err.println("Waiting for game tree to become available.");
            Thread.sleep(1000);
        }

        // Main game loop
        while (true) {

            // Get the strategy request for the current game state
            PlayStrategy myStrategy = this.getStrategyRequest();
            if (myStrategy == null) // Game was terminated by an outside event
                break;

            boolean playComplete = false; // Flag to check if the strategy was valid

            // Loop until a valid strategy is found
            while (!playComplete) {

                // Printing Separator
                System.out.println("*******************************************************");

                System.out.println("TournamentMixedNEStrategy");

                // Initialize variables
                GameNode rootNode = tree.getRootNode();

                // Get the number of decisions for each player
                int numberOfDecisionsP1 = rootNode.numberOfChildren();
                Iterator<GameNode> childrenNodesP1 = rootNode.getChildren();
                int numberOfDecisionsP2 = childrenNodesP1.next().numberOfChildren();

                // Initialize arrays to hold actions and utilities for both players
                String[] actionsP1 = new String[numberOfDecisionsP1];
                String[] actionsP2 = new String[numberOfDecisionsP2];
                int[][] utility1 = new int[numberOfDecisionsP1][numberOfDecisionsP2];
                int[][] utility2 = new int[numberOfDecisionsP1][numberOfDecisionsP2];
                GameNode childNodeP1;
                GameNode childNode2;
                // Initialize the index for the actions of Player 1 and 2
                int j;
                int i = 0;

                // Find payoffs for each possible pure strategy profile
                System.out.println("Assuming normal form game, see comment in the code.");
                // Get an iterator over the children of the root node
                childrenNodesP1 = rootNode.getChildren();

                // Iterate over each child node of the root node
                while (childrenNodesP1.hasNext()) {

                    // Get the next child node
                    childNodeP1 = childrenNodesP1.next();

                    // Store the label of the child node as an action for Player 1
                    actionsP1[i] = childNodeP1.getLabel();

                    // Reset the counter for Player 2's actions
                    j = 0;

                    // Get an iterator over the children of the current child node
                    Iterator<GameNode> childrenNodesP2 = childNodeP1.getChildren();

                    // Iterate over each action of Player 2
                    while (childrenNodesP2.hasNext()) {

                        // Get the next child node
                        childNode2 = childrenNodesP2.next();

                        // If this is the first iteration, store the label of the child node as an action for Player 2
                        if (i == 0)
                            actionsP2[j] = childNode2.getLabel(); // Assuming normal form game (perfect information not supported)

                        // Store the payoff for Player 1 for the current strategy profile
                        utility1[i][j] = childNode2.getPayoffP1();

                        // Store the payoff for Player 2 for the current strategy profile
                        utility2[i][j] = childNode2.getPayoffP2();

                        // Increment the counter for Player 2's actions
                        j++;
                    }

                    // Increment the counter for Player 1's actions
                    i++;
                }

                // Compute Mixed Nash Equilibrium
                double[] strategyP1 = computeMixedStrategy(utility2);
                double[] strategyP2 = computeMixedStrategy(utility1);

                // Print the strategies for both players
                showStrategy(1, strategyP1, actionsP1);
                showStrategy(2, strategyP2, actionsP2);

                // Assuming again it's n x m
                for (i = 0; i < actionsP1.length; i++) {
                    myStrategy.put(actionsP1[i], strategyP1[i]);
                }
                for (j = 0; j < actionsP2.length; j++) {
                    myStrategy.put(actionsP2[j], strategyP2[j]);
                }

                // Try to provide the computed strategy for the game
                try {
                    this.provideStrategy(myStrategy);
                    playComplete = true;
                } catch (InvalidStrategyException e) {
                    System.err.println("Invalid strategy: " + e.getMessage());
                    e.printStackTrace(System.err);
                }

            }

            // Try to provide the computed strategy for the game
            try {
                this.provideStrategy(myStrategy);
                playComplete = true;
            } catch (InvalidStrategyException e) {
                System.err.println("Invalid strategy: " + e.getMessage());
                e.printStackTrace(System.err);
            }
        }
    }

    private double[] computeMixedStrategy(int[][] utilities) {
        int n = utilities.length; // Number of strategies for Player 1
        int m = utilities[0].length; // Number of strategies for Player 2

        // Coefficients for the objective function
        double[] c = new double[n];
        Arrays.fill(c, 1.0);

        // Constraints: Expected utility for each strategy of Player 2 must be equal
        double[][] A = new double[m + 1][n];
        double[] b = new double[m + 1];

        for (int j = 0; j < m; j++) {
            for (int i = 0; i < n; i++) {
                A[j][i] = utilities[i][j];
            }
            b[j] = 0;
        }

        // Sum of probabilities must be equal to 1
        for (int i = 0; i < n; i++) {
            A[m][i] = 1.0;
        }
        b[m] = 1.0;

        // Lower bounds for the probabilities
        double[] lb = new double[n];
        for (int i = 0; i < n; i++) {
            lb[i] = 0.0;
        }

        // Set up and solve the linear programming problem
        LinearProgram lp = new LinearProgram(c);
        lp.setMinProblem(false);
        for (int i = 0; i < b.length; i++) {
            lp.addConstraint(new LinearSmallerThanEqualsConstraint(A[i], b[i], "c" + i));
        }
        lp.setLowerbound(lb);

        double[] solution = solveLP(lp);
        if (solution == null) {
            throw new RuntimeException("No solution found for the linear program");
        }

        return solution;
    }

    public double[] solveLP(LinearProgram lp) {
        LinearProgramSolver solver = SolverFactory.newDefault();
        double[] x = solver.solve(lp);
        return x;
    }

    // Method to extract the label from a string
    public String showLabel(String label) {
        return label.substring(label.lastIndexOf(':') + 1);
    }

    // Method to print the strategy for a player
    public void showStrategy(int P, double[] strategy, String[] labels) {
        System.out.println("Strategy Player " + P + ":");
        for (int i = 0; i < labels.length; i++) {
            System.out.println("   " + strategy[i] + ":" + showLabel(labels[i]));
        }
    }
}
