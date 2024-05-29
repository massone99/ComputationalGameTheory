package play;

import gametree.GameNode;
import scpsolver.constraints.LinearBiggerThanEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;

import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class Quiz1_3Strategy extends Strategy {

    static LinearProgram lp; // The linear program
    static double[] x; // The solution to the linear program
    static List<Integer> activeP1; // The indices of the remaining actions for player 1
    static List<Integer> activeP2; // The indices of the remaining actions for player 2

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

                System.out.println("0");

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

                System.out.println("ok");

                // Find payoffs for each possible pure strategy profile
                System.out.println("Assuming normal form game, see comment in the code.");
                Iterator<GameNode> childrenNodes1 = rootNode.getChildren();

                System.out.println("okok");

                // Compute U1 and U2
                while (childrenNodes1.hasNext()) {
                    childNode1 = childrenNodes1.next();
                    labelsP1[i] = childNode1.getLabel();
                    j = 0;
                    Iterator<GameNode> childrenNodes2 = childNode1.getChildren();
                    while (childrenNodes2.hasNext()) {
                        childNode2 = childrenNodes2.next();
                        if (i == 0)
                            labelsP2[j] = childNode2.getLabel(); // Assuming normal form game (perfect information not supported)
                        U1[i][j] = childNode2.getPayoffP1();
                        U2[i][j] = childNode2.getPayoffP2();
                        j++;
                    }
                    i++;
                }

                System.out.println("okokok");

                // Reduce
                if( n1 != 2 || n2 != 2 ) {
                    setLPQuiz1(U1, U2, labelsP1, labelsP2);
                }

                // If 2x2 compute nash
                if( n1 == 2 && n2 == 2 ) {
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
                }

                // Now actually play a random strategy

                System.out.println("okokokok");

                // Get iterator of number of moves per turn
                Iterator<Integer> iterator = tree.getValidationSet().iterator();
                // Get iterator of moves already in the strategy
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
                            return;
                        }
                        myStrategy.put(keys.next(), moves[i]);
                    }
                }

            }
        }
    }


    private static boolean isDominated(int[][] payoffs, List<Integer> activeSelf, List<Integer> activeOpponent, int index, boolean isPlayerOne) {

        if (activeSelf.size() == 1){
            return false;
        }

        int minValue = 0;

        for (int i = 0 ; i < payoffs.length ; i++){
            if (!activeSelf.contains(i) || i == index){
                continue;
            }

            for (int j = 0 ; j < payoffs[i].length ; j++){

                if (!activeOpponent.contains(j)){
                    continue;
                }

                int newMin = payoffs[i][j];

                if (newMin > 0){
                    newMin = 0;
                }

                minValue = minValue > newMin ? newMin : minValue;
            }
        }

        if (minValue < 0){
            minValue = -1 * minValue;
        }

        System.out.println("Minimum value for " + (isPlayerOne ? "row" : "column") + " " + index + ": " + minValue);
        double[] c = new double[activeSelf.size() - 1];
        Arrays.fill(c, 1); // Objective: minimize the sum of probabilities

        double[] lb = new double[activeSelf.size() - 1];
        Arrays.fill(lb, 0);

        lp = new LinearProgram(c);
        lp.setMinProblem(true);

        for (int j = 0; j < activeOpponent.size(); j++) {

            Boolean skippedIndex = false;
            double[] constraintCoefficients = new double[activeSelf.size() - 1];
            for (int i = 0; i < activeSelf.size(); i++) {
                if (i == index){
                    skippedIndex = true;
                    continue;
                }

                if (!skippedIndex){
                    constraintCoefficients[i] = payoffs[activeSelf.get(i)][activeOpponent.get(j)] + minValue;
                }
                else{
                    constraintCoefficients[i-1] = payoffs[activeSelf.get(i)][activeOpponent.get(j)] + minValue;
                }
            }

            System.out.println("Adding constraint for " + (isPlayerOne ? "row" : "column") + " " + index);
            lp.addConstraint(new LinearBiggerThanEqualsConstraint(constraintCoefficients, payoffs[activeSelf.get(index)][activeOpponent.get(j)] + minValue, "c" + j));
        }

        lp.setLowerbound(lb);

        boolean solvable = solveLP();

        return solvable && lp.evaluate(x) < 1 && lp.evaluate(x) >= 0;
    }

    public static boolean solveLP() {
        LinearProgramSolver solver  = SolverFactory.newDefault();
        x = solver.solve(lp);
        return x != null;
    }

    public static void setLPQuiz1(int[][] player1Payoffs, int[][] player2Payoffs, String[] p1Labels, String[] p2Labels) {
        activeP1 = new ArrayList<>();
        activeP2 = new ArrayList<>();

        for (int i = 0; i < p1Labels.length; i++) {
            activeP1.add(i);
        }

        List<Integer> activeP2 = new ArrayList<>();
        for (int j = 0; j < p2Labels.length; j++) {
            activeP2.add(j);
        }

        boolean change;
        do {
            change = false;


            for (int i = 0; i < activeP1.size(); i++) {
                if (isDominated(player1Payoffs, activeP1, activeP2, i, true)) {
                    activeP1.remove(i);
                    change = true;
                    break;
                }
            }

            if (change) continue;

            for (int i = 0; i < activeP2.size(); i++) {
                if (isDominated(transpose(player2Payoffs), activeP2, activeP1, i, false)) {
                    activeP2.remove(i);
                    change = true;
                    break;
                }
            }
        } while (change);

        System.out.println("Remaining Player 1 actions: " + indicesToLabels(activeP1, p1Labels));
        System.out.println("Remaining Player 2 actions: " + indicesToLabels(activeP2, p2Labels));
    }

    private static int[][] transpose(int[][] matrix) {
        int[][] transposed = new int[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                transposed[j][i] = matrix[i][j];
            }
        }
        return transposed;
    }

    private static List<String> indicesToLabels(List<Integer> indices, String[] labels) {
        List<String> result = new ArrayList<>();
        for (int index : indices) {
            result.add(labels[index]);
        }
        return result;
    }


}
