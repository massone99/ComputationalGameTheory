package play;

import gametree.GameNode;
import play.exception.InvalidStrategyException;

import java.util.Iterator;

public class MixedNEStrategy extends Strategy {

    @Override
    public void execute() throws InterruptedException {

        // Wait until the game tree is available
        while( !this.isTreeKnown() ) {
            System.err.println("Waiting for game tree to become available.");
            Thread.sleep(1000);
        }

        // Main game loop
        while( true ) {

            // Get the strategy request for the current game state
            PlayStrategy myStrategy = this.getStrategyRequest();
            if( myStrategy == null ) // Game was terminated by an outside event
                break;

            boolean playComplete = false; // Flag to check if the strategy was valid

            // Loop until a valid strategy is found
            while( !playComplete ) {

                // Printing Separator
                System.out.println("*******************************************************");

                // Initialize variables
                GameNode rootNode = tree.getRootNode();

                // Get the number of decisions for each player
                int numberOfDecisionsP1 = rootNode.numberOfChildren();
                int numberOfDecisionsP2 = rootNode.getChildren().next().numberOfChildren();

                // Check if the game is 2x2
                if( numberOfDecisionsP1 != 2 || numberOfDecisionsP2 != 2 ) {
                    System.err.println("PANIC: Game is not 2x2.");
                    return;
                }

                // Initialize arrays to hold actions and utilities for both players
                String[] actionsP1 = new String[numberOfDecisionsP1];
                String[] actionsP2 = new String[numberOfDecisionsP2];
                int[][] utility1 = new int[numberOfDecisionsP1][numberOfDecisionsP2];
                int[][] utility2 = new int[numberOfDecisionsP1][numberOfDecisionsP2];
                GameNode childNodeP1, childNode2;
                // Initialize the index for the actions of Player 1 and 2
                int j, i = 0;

                // Find payoffs for each possible pure strategy profile
                System.out.println("Assuming normal form game, see comment in the code.");
                // Get an iterator over the children of the root node
                Iterator<GameNode> childrenNodesP1 = rootNode.getChildren();

                // Iterate over each child node of the root node
                while( childrenNodesP1.hasNext() ) {

                    // Get the next child node
                    childNodeP1 = childrenNodesP1.next();

                    // Store the label of the child node as an action for Player 1
                    actionsP1[i] = childNodeP1.getLabel();

                    // Reset the counter for Player 2's actions
                    j = 0;

                    // Get an iterator over the children of the current child node
                    Iterator<GameNode> childrenNodesP2 = childNodeP1.getChildren();

                    // Iterate over each action of Player 2
                    while( childrenNodesP2.hasNext() ) {

                        // Get the next child node
                        childNode2 = childrenNodesP2.next();

                        // If this is the first iteration, store the label of the child node as an action for Player 2
                        if ( i==0 )
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
                double[] strategyP1 = new double[numberOfDecisionsP1];
                double[] strategyP2 = new double[numberOfDecisionsP2];

                // Compute strategy for Player 1
                if ((utility2[0][0] - utility2[1][0] + utility2[1][1] - utility2[0][1]) == 0) {
                    System.out.println("No Nash equilibrium for player 1.");
                    // Random because either any strategy is a Nash equilibrium
                    // or none is a Nash equilibrium
                    strategyP1[0] = 0.5;
                    strategyP1[1] = 0.5;
                } else {
                    strategyP1[0] = (double) (utility2[1][1] - utility2[1][0]) / (utility2[0][0] - utility2[1][0] + utility2[1][1] - utility2[0][1]);
                    System.out.println("Strategy 1: " + strategyP1[0]);
                    strategyP1[1] = 1 - strategyP1[0];
                }

                // Compute strategy for Player 2
                if ((utility1[0][0] - utility1[1][0] + utility1[1][1] - utility1[0][1]) == 0) {
                    System.out.println("No Nash equilibrium for player 2.");
                    // Random because either any strategy is a Nash equilibrium
                    // or none is a Nash equilibrium
                    strategyP2[0] = 0.5;
                    strategyP2[1] = 0.5;
                } else {
                    strategyP2[0] = (double) (utility1[1][1] - utility1[0][1]) / (utility1[0][0] - utility1[0][1] + utility1[1][1] - utility1[1][0]);
                    System.out.println("Strategy 2: " + strategyP1[0]);
                    strategyP2[1] = 1 - strategyP2[0];
                }

                // Print the strategies for both players
                showStrategy(1,strategyP1,actionsP1);
                showStrategy(2,strategyP2,actionsP2);

                // Assuming again it's 2x2
                for ( i = 0; i < actionsP1.length; i++ ) {
                    myStrategy.put(actionsP1[i], strategyP1[i]);
                    myStrategy.put(actionsP2[i], strategyP2[i]);
                }

                // Try to provide the computed strategy for the game
                try{
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
                System.err.println("Invalid strategy: " + e.getMessage());;
                e.printStackTrace(System.err);
            }
        }
    }

    // Method to extract the label from a string
    public String showLabel(String label) {
        return label.substring(label.lastIndexOf(':')+1);
    }

    // Method to print the strategy for a player
    public void showStrategy(int P, double[] strategy, String[] labels) {
        System.out.println("Strategy Player " + P + ":");
        for (int i = 0; i<labels.length; i++) System.out.println("   " + strategy[i] + ":" + showLabel(labels[i]));
    }

}