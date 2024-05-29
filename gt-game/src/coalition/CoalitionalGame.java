package coalition;

import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import scpsolver.constraints.LinearBiggerThanEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class CoalitionalGame {
    private static final int BASE = 2;
    private static final double ROUNDING_SCALE = 100.0;
    private static final int ASCII_OFFSET = 64;
    public final double[] payoffs;
    public final int nPlayers;
    public final Map<String, Integer> idMap; // New map to store the mapping from each id to an integer
    public String[] ids;
    private double quota;

    public CoalitionalGame(double[] payoffs) {
        this.payoffs = payoffs;
        this.nPlayers = (int) (Math.log(payoffs.length) / Math.log(2));
        setPlayersID();
        idMap = new HashMap<>(); // Initialize the map
        for (int i = 0; i < nPlayers; i++) {
            idMap.put(ids[i], i); // Map each id to an integer
        }
    }

    public static void main(String[] args) {
        // double[] v1 = {0.0, 0.0, 3.0, 8.0, 4.0, 5.0, 7.0, 10.0};
        double[] v1 = readVectorFromFile("gt-game/resources/EC1.txt");
        CoalitionalGame c = new CoalitionalGame(v1);
        c.showGame();
        for (int j = 0; j < c.nPlayers; j++) {
            System.out.println("*********** Permutations without player " + c.ids[c.nPlayers - 1 - j] + " ***********");
            for (int i = 0; i < c.nPlayers; i++) {
                System.out.print("With " + i + " players: ");
                c.permutation(0, i, j, 0);
                System.out.println();
            }
        }
        Map<String, Double> shapleyValues = c.computeShapleyValues();
        boolean isShapleyInCore = c.isShapleyInCore(new ArrayList<>(shapleyValues.values()),
                new ArrayList<>(c.idMap.values()), c.nPlayers);
        System.out.println("Is Shapley value in the core? " + isShapleyInCore);
        List<Long> allCoalitions = c.getAllCoalitions();
        List<Set<String>> allCoalitionsLabels = new ArrayList<>();

        for (Long coalition : allCoalitions) {
            // Convert the coalition to binary format
            String binaryCoalition = Long.toBinaryString(coalition);

            // Get the labels of the players in the coalition
            ArrayList<Integer> playersInCoalition = c.getSet(coalition);
            List<String> labels = c.getAgentLabels(playersInCoalition);

            // Get the payoff of the coalition
            double payoff = c.getPayoffOfCoalition(coalition);

            // Print the coalition in binary format, the labels, and the payoff
            //System.out.println("Coalition: " + binaryCoalition + ", Labels: " + labels + ", Payoff: " + payoff);
        }

        //System.out.println("**************COALITIONS PAYOFFS**************");
        double[][] coalitionPayoffs = c.getAllCoalitionPayoffs();
        for (int i = 0; i < coalitionPayoffs.length; i++) {
            //System.out.println("Coalition: " + c.returnShowSetString(allCoalitions.get(i)) + ", Payoffs: " + Arrays.toString(coalitionPayoffs[i]));
        }

        System.out.println("**************CORE**************");
        List<double[]> core = c.findCore();
        for (double[] x : core) {
            System.out.println(Arrays.toString(x));
        }
    }

    private static double[] readVectorFromFile(String filename) {
        List<Double> list = new ArrayList<>();
        try (Scanner scanner = new Scanner(new File(filename))) {
            while (scanner.hasNextLine()) {
                list.add(Double.parseDouble(scanner.nextLine()));
            }
        } catch (FileNotFoundException ignored) {
        }
        double[] array = new double[list.size()];
        for (int i = 0; i < list.size(); i++) {
            array[i] = list.get(i);
        }
        return array;
    }

    /**
     * Compute the binomial coefficient C(n, k) = n! / (k! * (n - k)!)
     *
     * @param n
     * @param k
     * @return the binomial coefficient
     */
    private static long binomialCoefficient(int n, int k) {
        long res = 1;
        if (k > n - k) {
            k = n - k;
        }
        for (int i = 0; i < k; ++i) {
            res *= (n - i);
            res /= (i + 1);
        }
        return res;
    }

    public List<double[]> findCore() {
        List<double[]> core = new ArrayList<>();
        int n = nPlayers;

        // Iterate over all possible payoff vectors
        for (double[] x : getAllCoalitionPayoffs()) {
            boolean isInCore = true;

            // Check the core condition for each coalition
            for (Long coalition : getAllCoalitions()) {
                double coalitionValue = getValueOfCoalition(longToSet(coalition));
                double payoffSum = 0;

                // Iterate over each player
                for (int player = 0; player < n; player++) {
                    // Check if the player is in the coalition
                    if ((coalition & (1L << player)) != 0) {
                        // If the player is in the coalition, add their payoff to the total
                        payoffSum += x[player];
                    }
                }

                if (payoffSum < coalitionValue) {
                    isInCore = false;
                    break;
                }
            }

            if (isInCore) {
                core.add(x);
            }
        }

        return core;
    }

    public Set<Integer> longToSet(long coalition) {
        Set<Integer> set = new HashSet<>();
        int playerIndex = 0;
        while (coalition > 0) {
            if ((coalition & 1) == 1) {
                set.add(playerIndex);
            }
            coalition >>= 1;
            playerIndex++;
        }
        return set;
    }

    public double[][] getAllCoalitionPayoffs() {
        // Get all possible coalitions
        List<Long> allCoalitions = getAllCoalitions();

        // Initialize the 2D array
        double[][] coalitionPayoffs = new double[allCoalitions.size()][nPlayers];

        // Iterate over all coalitions
        for (int i = 0; i < allCoalitions.size(); i++) {
            Long coalition = allCoalitions.get(i);

            // Get the players in the coalition
            ArrayList<Integer> playersInCoalition = getSet(coalition);

            // Get the payoff of the coalition
            double payoff = getPayoffOfCoalition(coalition);

            // Iterate over all players
            for (int j = 0; j < nPlayers; j++) {
                // If the player is part of the coalition, store the payoff in the array
                if (playersInCoalition.contains(j)) {
                    coalitionPayoffs[i][j] = payoff;
                }
            }
        }

        // Initialize the inverted array
        double[][] invertedCoalitionPayoffs = new double[allCoalitions.size()][nPlayers];

        // Iterate over all rows
        for (int i = 0; i < coalitionPayoffs.length; i++) {
            // Iterate over all columns in reverse order
            for (int j = 0; j < coalitionPayoffs[i].length; j++) {
                // Fill the inverted array with the values from the original array in reverse order
                invertedCoalitionPayoffs[i][j] = coalitionPayoffs[i][coalitionPayoffs[i].length - j - 1];
            }
        }

        return invertedCoalitionPayoffs;
    }

    public double getPayoffOfCoalition(long coalition) {
        // Convert the coalition to an integer
        int index = (int) coalition;

        // Return the payoff of the coalition
        return payoffs[index];
    }

    /**
     * Generate all possible coalitions in the game.
     *
     * @return a list of all possible coalitions (represented as sets of player
     * indices)
     */
    public List<Long> getAllCoalitions() {
        List<Long> allCoalitions = new ArrayList<>();
        long totalCoalitions = (long) Math.pow(2, nPlayers) - 1; // Total possible coalitions excluding the empty set
        for (long i = 1; i <= totalCoalitions; i++) {
            allCoalitions.add(i);
        }
        return allCoalitions;
    }

//    public List<double[]> calculateCore() {
//        List<double[]> core = new ArrayList<>();
//        int n = nPlayers;
//
//        // Iterate over all possible payoff vectors
//        for (double[] x : getAllPayoffVectors(n)) {
//            boolean isInCore = true;
//
//            // Check the core condition for each coalition
//            for (Set<Integer> coalition : getAllCoalitions()) {
//                double coalitionValue = getValueOfCoalition(coalition);
//                double payoffSum = 0;
//
//                for (int player : coalition) {
//                    payoffSum += x[player]; // Assuming x is the payoff vector
//                }
//
//                if (payoffSum < coalitionValue) {
//                    isInCore = false;
//                    break;
//                }
//            }
//
//            if (isInCore) {
//                core.add(x);
//            }
//        }
//
//        return core;
//    }

    public List<String> getAgentLabels(List<Integer> agentIndices) {
        List<String> agentLabels = new ArrayList<>();
        for (Integer index : agentIndices) {
            agentLabels.add(ids[index]);
        }
        return agentLabels;
    }

    public boolean isShapleyInCore(List<Double> shapleyValues, List<Integer> ids, int nPlayers) {
        List<Integer> players = new ArrayList<>();
        for (int i = 0; i < nPlayers; i++) {
            players.add(i);
        }
        Collections.reverse(ids);
        return checkAllSubsets(players, shapleyValues, ids, new ArrayList<>(), 0);
    }

    private boolean checkAllSubsets(List<Integer> players, List<Double> shapleyValues, List<Integer> ids,
                                    List<Integer> current, int index) {
        if (index == players.size()) {
            double coalitionShapleyValue = 0.0;
            // Calculate the total Shapley value for the current coalition
            for (Integer player : current) {
                coalitionShapleyValue += shapleyValues.get(ids.get(player));
            }
            // Check if the total Shapley value is less than the value of the coalition
            // If it is, return false (the Shapley value condition is not satisfied)
            return !(coalitionShapleyValue < getValueOfCoalition(new HashSet<>(current)));
        } else {
            // Add the current player to the coalition
            current.add(players.get(index));
            // Recursively check the next player
            if (!checkAllSubsets(players, shapleyValues, ids, current, index + 1)) {
                return false;
            }
            // Remove the current player from the coalition
            current.remove(current.size() - 1);
            // Recursively check the next player without the current player in the coalition
            return checkAllSubsets(players, shapleyValues, ids, current, index + 1);
        }
        // If all subsets satisfy the Shapley value condition, return true
    }

    public double getValueOfCoalition(Set<Integer> coalition) {
        int index = 0;
        for (Integer player : coalition) {
            index |= 1 << player;
        }
        return payoffs[index];
    }

    public Map<String, Double> computeShapleyValues() {
        System.out.println("*********** Computing Shapley Values ***********");
        Map<String, Double> shapleyValues = new TreeMap<>(); // TreeMap will sort the keys in natural order

        for (int i = 0; i < nPlayers; i++) {
            double shapleyValue = computeShapleyValueForPlayer(i);
            shapleyValues.put(ids[i], shapleyValue);
        }

        printShapleyValues(shapleyValues);

        return shapleyValues;
    }

    private double computeShapleyValueForPlayer(int playerIndex) {
        double shapleyValue = 0;

        for (int j = 0; j < Math.pow(2, nPlayers); j++) {
            if ((j & (1 << playerIndex)) != 0) {
                shapleyValue += computeGainForPlayerInCoalition(playerIndex, j);
            }
        }

        return Math.round(shapleyValue * ROUNDING_SCALE) / ROUNDING_SCALE; // Round to the second decimal
    }

    private double computeGainForPlayerInCoalition(int playerIndex, int coalition) {
        // Compute the coalition without player i
        int coalitionWithoutPlayer = coalition - (1 << playerIndex);
        // Compute the gain of player i when joining the coalition
        double gain = payoffs[coalition] - payoffs[coalitionWithoutPlayer];
        // Compute the weight of the coalition
        double weight = 1.0 / (nPlayers * binomialCoefficient(nPlayers - 1, Integer.bitCount(coalition) - 1));
        // Update the Shapley value
        double weightedGain = gain * weight;

        System.out.println(returnShowSetString(coalitionWithoutPlayer) + " (" + payoffs[coalitionWithoutPlayer] + ") -> "
                + returnShowSetString(coalition)
                + " (" + payoffs[coalition] + ") ***** gain = " + gain + " (" + weight + ")");

        return weightedGain;
    }

    private void printShapleyValues(Map<String, Double> shapleyValues) {
        // Print Shapley values in alphabetical order based on the ids
        for (Map.Entry<String, Double> entry : shapleyValues.entrySet()) {
            System.out.println("Shapley Value for " + entry.getKey() + " : " + entry.getValue());
        }
    }

    public void setPlayersID() {
        int c = ASCII_OFFSET;
        ids = new String[nPlayers];
        for (int i = nPlayers - 1; i >= 0; i--) {
            c++;
            ids[i] = (String.valueOf((char) c));
        }
    }

    public void showGame() {
        System.out.println("*********** Coalitional Game ***********");
        for (int i = 0; i < payoffs.length; i++) {
            showSet(i);
            System.out.println(" (" + payoffs[i] + ")");
        }
    }

    public String returnShowSetString(long v) {
        StringBuilder sb = new StringBuilder();
        boolean showPlayerID = true;
        int power;
        sb.append("{");
        int cnt = 0;
        for (int i = 0; i < nPlayers; i++) {
            power = nPlayers - (i + 1);
            if (showPlayerID) {
                if (inSet(i, v)) {
                    if (cnt > 0)
                        sb.append(",");
                    cnt++;
                    sb.append(ids[power]);
                }
            } else {
                if (cnt > 0)
                    sb.append(",");
                cnt++;
                if (inSet(i, v))
                    sb.append(1);
                else
                    sb.append(0);
            }
        }
        sb.append("}");
        return sb.toString();
    }

    public void showSet(long v) {
        int power;
        StringBuilder output = new StringBuilder("{");

        for (int i = 0; i < nPlayers; i++) {
            power = nPlayers - (i + 1);
            if (inSet(i, v)) {
                if (output.length() > 1) {
                    output.append(",");
                }
                output.append(ids[power]);
            }
        }

        output.append("}");
        System.out.print(output.toString());
    }

    /**
     * This method checks if a player is in a coalition.
     *
     * @param i The index of the player. This is a zero-based index where 0 represents the first player.
     * @param v The coalition represented as a long. Each bit in the binary representation of v represents a player, where a 1 means the player is in the coalition and a 0 means the player is not.
     * @return A boolean value indicating whether the player is in the coalition. It returns true if the player is in the coalition and false otherwise.
     * <p>
     * The method works by first calculating the power of 2 corresponding to the player's position. This is done by subtracting the player's index (plus one) from the total number of players. This gives us the position of the player in the binary representation of the coalition (from right to left, starting from 0).
     * <p>
     * The method then checks if the binary representation of the coalition includes a 1 at the position corresponding to the player. This is done by right shifting the bits of v by the calculated power and then performing a bitwise AND operation with 1. The right shift operation effectively moves the bit at the player's position to the least significant bit position. The bitwise AND operation then checks if this bit is 1.
     * <p>
     * If the bit at the player's position is 1, the method returns true, indicating that the player is in the coalition. If the bit is 0, the method returns false.
     */
    public boolean inSet(int i, long v) {
        int power = nPlayers - (i + 1);
        return ((v >> power) & 1) == 1;
    }

    public ArrayList<Integer> getSet(long v) {
        ArrayList<Integer> players = new ArrayList<>();
        int power;
        long vi;
        long div;
        long mod;
        for (int i = 0; i < nPlayers; i++) {
            power = nPlayers - (i + 1);
            vi = (long) Math.pow(2, power);
            div = v / vi;
            mod = div % 2;
            if (mod == 1)
                players.add(power);
        }
        return players;
    }

    public void permutation(int currentPlayerIndex, int remainingPlayers, int excludedPlayerIndex,
                            long currentSetValue) {
        // Initialize the accumulated value to 0
        long accumulatedValue = 0;

        // If there are no remaining players, display the current set value
        if (remainingPlayers == 0) {
            showSet(currentSetValue);
        } else {
            // Initialize the offset to 0
            int offset = 0;

            // Calculate the offset based on the excluded player index and current player
            // index
            if (excludedPlayerIndex < currentPlayerIndex)
                offset = nPlayers - currentPlayerIndex;
            else
                offset = nPlayers - currentPlayerIndex - 1;

            // If the offset equals the remaining players
            if (offset == remainingPlayers) {
                // Loop through the players from the current player index
                for (int i = currentPlayerIndex; i < nPlayers; i++) {
                    // If the current player is not the excluded player
                    if (i != excludedPlayerIndex)
                        // Add the value of the current player to the accumulated value
                        accumulatedValue += (long) Math.pow(2, nPlayers - (i + 1));
                }

                // Add the accumulated value to the current set value
                currentSetValue = currentSetValue + accumulatedValue;

                // Display the current set value
                showSet(currentSetValue);
            } else {
                // If the current player is not the excluded player, generate the next
                // permutation
                if (currentPlayerIndex != excludedPlayerIndex)
                    permutation(currentPlayerIndex + 1, remainingPlayers - 1, excludedPlayerIndex,
                            currentSetValue + (long) Math.pow(2, nPlayers - (currentPlayerIndex + 1)));

                // Generate the next permutation
                permutation(currentPlayerIndex + 1, remainingPlayers, excludedPlayerIndex, currentSetValue);
            }
        }
    }

}
