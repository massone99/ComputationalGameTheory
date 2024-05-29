package play.quiz2;

import gametree.GameNode;
import play.LinearProgramming;
import play.NormalFormGame;
import play.PlayStrategy;
import play.Strategy;
import play.exception.InvalidStrategyException;
import scpsolver.constraints.LinearEqualsConstraint;
import scpsolver.constraints.LinearSmallerThanEqualsConstraint;
import scpsolver.problems.LinearProgram;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class Quiz2_1Strategy extends Strategy {
    private static int count = 0;

    // Method to transpose a given matrix
    public static int[][] transpose(int[][] m) {
        // Create a new matrix with dimensions swapped for transposition
        int[][] transposed = new int[m[0].length][m.length];
        // Iterate through each element of the original matrix
        for ( int i = 0; i < m.length; i++ ) {
            for ( int j = 0; j < m[0].length; j++ )
                // Swap the row and column indices to transpose the element
                transposed[j][i] = m[i][j];
        }
        // Return the transposed matrix
        return transposed;
    }


    // Method to generate all possible subsets of size 'k' from a set of size 'n'
    // 'j' is the current index, 'k' is the size of the subset to be generated,
    // 'n' is the size of the original set, and 'p' is a boolean array indicating
    // which elements are included in the original set
    public static List<boolean[]> getSubSets(int j, int k, int n, boolean[] incl) {
        boolean[] sub = new boolean[n]; // Create a boolean array to represent subsets
        List<boolean[]> subset = new ArrayList<>(); // Initialize a list to store subsets

        // Base case: if k is 0, add an empty subset to the list
        if ( k == 0 ) {
            subset.add(sub);
        } else {
            int c = 0; // Count the number of true values in p
            for ( int i = j; i < incl.length; i++ ) {
                if ( incl[i] )
                    c++;
            }

            // If op equals k, create a subset containing all true values in p
            if ( c == k ) {
                for ( int i = 0; i < sub.length; i++ ) {
                    if ( incl[i] )
                        sub[i] = true;
                }
                subset.add(sub);
            } else {
                // Recursively generate subsets with and without including the current element
                if ( incl[j] ) {
                    List<boolean[]> s1 = getSubSets(j + 1, k - 1, n, incl);
                    for (boolean[] booleans : s1) {
                        sub = booleans;
                        sub[j] = true;
                        subset.add(sub);
                    }
                }
                List<boolean[]> s0 = getSubSets(j + 1, k, n, incl);
                for (boolean[] booleans : s0) {
                    sub = booleans;
                    sub[j] = false;
                    subset.add(sub);
                }
            }
        }
        return subset; // Return the list of subsets
    }

    // Main function to find Nash equilibrium with support
    public static void findNashAndPrint(int[][] payoffsP1, int[][] payoffsP2, boolean[] p1Sub, boolean[] p2Sub, String[] labelsP1, String[] labelsP2) {
        // Transpose player 2's payoffs matrix
        int[][] player2PayoffsTransposed = transpose(payoffsP2);

        // Get the number of strategies for both players
        int numStOfP1 = payoffsP1.length;
        int numStOfP2 = payoffsP1[0].length;

        // Calculate the size of the linear program
        int problemSize = numStOfP1 + numStOfP2 + 2;

        // Initialize arrays for the objective function coefficients, lower bounds, and constraint coefficients
        double[] c = new double[problemSize]; // Objective function coefficients
        Arrays.fill(c, 0.0); // Initialize with zeros

        double[] lb = new double[problemSize]; // Lower bounds
        Arrays.fill(lb, 0.0); // Initialize with zeros

        double[] BConstraint = new double[problemSize]; // Constraint coefficients
        Arrays.fill(BConstraint, 0.0); // Initialize with zeros

        // Initialize a new linear program
        LinearProgramming.lp = new LinearProgram(c);
        LinearProgramming.lp.setMinProblem(true); // Set it as a minimization problem

        // Loop through all strategies for both players
        for ( int k = 0; k < problemSize - 2; k++ ) {

            // Initialize constraint coefficients for each constraint
            double[] AConstraint = new double[problemSize];

            // Set up constraints based on whether the current strategy belongs to player 1 or player 2
            if ( k < numStOfP2 )
                AConstraint[problemSize - 2] = -1.0; // Constraint for player 1's strategy
            else
                AConstraint[problemSize - 1] = -1.0; // Constraint for player 2's strategy

            // Loop through all strategies again to set up payoff constraints
            for ( int i = 0; i < problemSize - 2; i++ ) {
                // Add constraints for player 1
                if ( i < numStOfP1 && k < numStOfP2 ) {
                    if ( p1Sub[i] )
                        AConstraint[i] = player2PayoffsTransposed[k][i];
                    else
                        AConstraint[i] = 0.0;
                }
                // Add constraints for player 2
                if ( i >= numStOfP1 && k >= numStOfP2 ) {
                    if ( p2Sub[i - numStOfP1] )
                        AConstraint[i] = payoffsP1[k - numStOfP2][i - numStOfP1];
                    else
                        AConstraint[i] = 0.0;
                }
            }

            // Determine whether to add an equality or inequality constraint based on strategy inclusion
            if ( k < numStOfP2 ) {
                if (p2Sub[k])
                    LinearProgramming.lp.addConstraint(new LinearEqualsConstraint(AConstraint, BConstraint[k], "c" + k));
                else
                    LinearProgramming.lp.addConstraint(new LinearSmallerThanEqualsConstraint(AConstraint, BConstraint[k], "c" + k));
            } else {
                if (p1Sub[k - numStOfP2])
                    LinearProgramming.lp.addConstraint(new LinearEqualsConstraint(AConstraint, BConstraint[k], "c" + k));
                else
                    LinearProgramming.lp.addConstraint(new LinearSmallerThanEqualsConstraint(AConstraint, BConstraint[k], "c" + k));
            }
        }

        // Set up equality constraints for the sum of probabilities for each player's strategies
        double[] eq1AConstraint = new double[problemSize];
        Arrays.fill(eq1AConstraint, 0.0);
        for ( int i = 0; i < p1Sub.length; i++ ) {
            if (p1Sub[i])
                eq1AConstraint[i] = 1.0;
        }

        double[] eq2AConstraint = new double[problemSize];
        Arrays.fill(eq2AConstraint, 0.0);
        for ( int i = 0; i < p2Sub.length; i++ ) {
            if ( p2Sub[i] )
                eq2AConstraint[i + numStOfP1] = 1.0;
        }

        // Add constraints for the sum of probabilities for player 1 and player 2
        int idx = problemSize - 2; // Index for player 1's sum constraint
        LinearProgramming.lp.addConstraint(new LinearEqualsConstraint(eq1AConstraint, 1.0, "c" + idx));

        idx = problemSize - 1; // Index for player 2's sum constraint
        LinearProgramming.lp.addConstraint(new LinearEqualsConstraint(eq2AConstraint, 1.0, "c" + idx));

        // Set lower bounds for the linear program
        LinearProgramming.lp.setLowerbound(lb);

        // Solve the linear program
        if ( LinearProgramming.solveLP() )
            // If the linear program is solved, print the results
            show(labelsP1, labelsP2);
    }


    private static void show(String[] player1Labels, String[] player2Labels) {
        count++;

        String fileName = "output.txt";

        try (PrintWriter out = new PrintWriter(new FileWriter(fileName, true))) {
            out.println("###############################################################");
            out.println("Player 1:");

            int length = player1Labels.length;
            for (int i = 0; i < player1Labels.length; i++)
                out.println(player1Labels[i] + " = " + String.format("%.2f", LinearProgramming.x[i]));

            out.println("###############################################################");
            out.println("Player 2:");

            for (int i = 0; i < player2Labels.length; i++)
                System.out.println(player2Labels[i] + " = " + String.format("%.2f", LinearProgramming.x[i + length]));

            out.println("The overall number of solutions is: " + count);
            out.println("###############################################################");
            out.println("");
        } catch (IOException e) {
            System.out.println("An error occurred w/ file: " + e.getMessage());
        }
    }


    public static void setLPQuiz2(int[][] player1Payoffs, int[][] player2Payoffs, String[] player1Labels, String[] player2Labels) {
        // Reset total solutions count
        count = 0;
        new java.io.File("output.txt").delete();

        // Determine the dimensions of the payoff matrices
        int maxSpace1 = player1Payoffs.length; // Number of rows
        int maxSpace2 = player1Payoffs[0].length; // Number of columns

        // Initialize boolean arrays to track subsets of strategies for each player
        boolean[] p1Sub = new boolean[maxSpace1]; // Player 1 subsets
        Arrays.fill(p1Sub, true); // Fill with true values
        boolean[] p2Sub = new boolean[maxSpace2]; // Player 2 subsets
        Arrays.fill(p2Sub, true); // Fill with true values

        // Loop over all possible support sizes
        for ( int supportSize = 1; supportSize <= (Math.min(maxSpace1, maxSpace2)); supportSize++) {
            // Generate all possible subsets for both players
            List<boolean[]> subsetsP1 = getSubSets(0, supportSize, maxSpace1, p1Sub); // Get subsets for player 1
            List<boolean[]> subsetsP2 = getSubSets(0, supportSize, maxSpace2, p2Sub); // Get subsets for player 2

            // Loop over all combinations of player 1 and player 2 subsets
            for (boolean[] supportP1 : subsetsP1) {
                for (boolean[] supportP2 : subsetsP2)
                    findNashAndPrint(player1Payoffs, player2Payoffs, supportP1, supportP2, player1Labels, player2Labels);
            }
        }
    }


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

            boolean playComplete = false;

            while( !playComplete ) {

                // Printing Separator
                System.out.println("*******************************************************");

                // Init variables
                GameNode rootNode = tree.getRootNode();

                int n1 = rootNode.numberOfChildren();
                int n2 = rootNode.getChildren().next().numberOfChildren();

                String[] labelsP1 = new String[n1];
                String[] labelsP2 = new String[n2];
                int[][] U1 = new int[n1][n2];
                int[][] U2 = new int[n1][n2];
                GameNode childNode1, childNode2;
                int j, i = 0;

                // Find payoffs for each possible pure strategy profile
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

                // Call the method to set up the linear program
                setLPQuiz2(U1, U2, labelsP1, labelsP2);

                // Get iterator of number of moves per turn
                Iterator<Integer> iterator = tree.getValidationSet().iterator();
                // Get iterator of strategies
                Iterator<String> keys = myStrategy.keyIterator();

                SecureRandom random = new SecureRandom();

                // Set random strategy
                while (iterator.hasNext()) {
                    double[] moves = new double[iterator.next()];
                    double sum = 0;

                    // For each move at this turn
                    // Compute a random probability such that they sum up to 1
                    for ( i = 0; i < moves.length - 1; i++ ) {
                        moves[i] = random.nextDouble();
                        while (sum + moves[i] >= 1)
                            moves[i] = random.nextDouble();
                        sum = sum + moves[i];
                    }
                    moves[moves.length - 1] = ((double) 1) - sum;

                    // Put in the strategy the probability for each action
                    for ( i = 0; i < moves.length; i++ ) {
                        if (!keys.hasNext()) {
                            System.err.println("PANIC: Strategy structure does not match the game.");
                            break;
                        }
                        myStrategy.put(keys.next(), moves[i]);
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

}
