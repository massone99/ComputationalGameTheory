package play;

import gametree.GameNode;
import play.exception.InvalidStrategyException;

import java.util.Iterator;

public class Quiz1_2Strategy extends Strategy {

    @Override
    public void execute() throws InterruptedException {

        while( !this.isTreeKnown() ) {
            System.err.println("Waiting for game tree to become available.");
            Thread.sleep(1000);
        }

        while( true ) {

            PlayStrategy myStrategy = this.getStrategyRequest();
            if( myStrategy == null ) // Game was terminated by an outside event
                break;

            boolean playComplete = false; // Strategy was not valid

            while( !playComplete ) {

                // Printing Separator
                System.out.println("*******************************************************");

                // Init variables
                GameNode rootNode = tree.getRootNode();

                int n1 = rootNode.numberOfChildren();
                int n2 = rootNode.getChildren().next().numberOfChildren();

                if( n1 != 2 || n2 != 2 ) {
                    System.err.println("PANIC: Game is not 2x2.");
                    return;
                }

                String[] labelsP1 = new String[n1];
                String[] labelsP2 = new String[n2];
                int[][] U1 = new int[n1][n2];
                int[][] U2 = new int[n1][n2];
                GameNode childNode1, childNode2;
                int j, i = 0;

                // Find payoffs for each possible pure strategy profile
                System.out.println("Assuming normal form game, see comment in the code.");
                Iterator<GameNode> childrenNodes1 = rootNode.getChildren();
                while( childrenNodes1.hasNext() ) {
                    childNode1 = childrenNodes1.next();
                    labelsP1[i] = childNode1.getLabel();
                    j = 0;
                    Iterator<GameNode> childrenNodes2 = childNode1.getChildren();
                    while( childrenNodes2.hasNext() ) {
                        childNode2 = childrenNodes2.next();
                        if ( i==0 )
                            labelsP2[j] = childNode2.getLabel(); // Assuming normal form game (perfect information not supported)
                        U1[i][j] = childNode2.getPayoffP1();
                        U2[i][j] = childNode2.getPayoffP2();
                        j++;
                    }
                    i++;
                }

                // Print
                new NormalFormGame(U1,U2,labelsP1,labelsP2).showGame();

                // Compute Mixed Nash Equilibrium
                double[] strategyP1 = new double[n1];
                double[] strategyP2 = new double[n2];

                // Player 1
                if ((U2[0][0] - U2[1][0] + U2[1][1] - U2[0][1]) == 0) {
                    System.out.println("No Nash equilibrium for player 1.");
                    // Random because either any strategy is a Nash equilibrium
                    // or none is a Nash equilibrium
                    strategyP1[0] = 0.5;
                    strategyP1[1] = 0.5;
                } else {
                    strategyP1[0] = (double) (U2[1][1] - U2[1][0]) / (U2[0][0] - U2[1][0] + U2[1][1] - U2[0][1]);
                    System.out.println("Strategy 1: " + strategyP1[0]);
                    strategyP1[1] = 1 - strategyP1[0];
                }

                // Player 2
                if ((U1[0][0] - U1[1][0] + U1[1][1] - U1[0][1]) == 0) {
                    System.out.println("No Nash equilibrium for player 2.");
                    // Random because either any strategy is a Nash equilibrium
                    // or none is a Nash equilibrium
                    strategyP2[0] = 0.5;
                    strategyP2[1] = 0.5;
                } else {
                    strategyP2[0] = (double) (U1[1][1] - U1[0][1]) / (U1[0][0] - U1[0][1] + U1[1][1] - U1[1][0]);
                    System.out.println("Strategy 2: " + strategyP1[0]);
                    strategyP2[1] = 1 - strategyP2[0];
                }

                showStrategy(1,strategyP1,labelsP1);
                showStrategy(2,strategyP2,labelsP2);

                // Assuming again it's 2x2
                for ( i = 0; i < labelsP1.length; i++ ) {
                    myStrategy.put(labelsP1[i], strategyP1[i]);
                    myStrategy.put(labelsP2[i], strategyP2[i]);
                }

                try{
                    this.provideStrategy(myStrategy);
                    playComplete = true;
                } catch (InvalidStrategyException e) {
                    System.err.println("Invalid strategy: " + e.getMessage());
                    e.printStackTrace(System.err);
                }

            }

            try {
                this.provideStrategy(myStrategy);
                playComplete = true;
            } catch (InvalidStrategyException e) {
                System.err.println("Invalid strategy: " + e.getMessage());;
                e.printStackTrace(System.err);
            }
        }
    }

    public String showLabel(String label) {
        return label.substring(label.lastIndexOf(':')+1);
    }

    public void showStrategy(int P, double[] strategy, String[] labels) {
        System.out.println("Strategy Player " + P + ":");
        for (int i = 0; i<labels.length; i++) System.out.println("   " + strategy[i] + ":" + showLabel(labels[i]));
    }

}
