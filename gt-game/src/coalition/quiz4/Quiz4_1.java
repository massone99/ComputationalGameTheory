package coalition.quiz4;

import scpsolver.constraints.Constraint;
import scpsolver.constraints.LinearBiggerThanEqualsConstraint;
import scpsolver.constraints.LinearEqualsConstraint;
import scpsolver.constraints.LinearSmallerThanEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.stream.Collectors;

public class Quiz4_1 {

    public final double[] v;
    public final int nPlayers;
    public final Map<String, Integer> idMap;
    public String[] ids;

    public Quiz4_1(double[] v) {
        this.v = v;
        this.nPlayers = (int) (Math.log(v.length) / Math.log(2));
        setPlayersID();
        idMap = new HashMap<>();
        for (int i = 0; i < nPlayers; i++) {
            idMap.put(ids[i], i);
        }
    }

    public static void main(String[] args) {
        double[] v1 = readVectorFromFile("gt-game/resources/C3.txt");
        Quiz4_1 c = new Quiz4_1(v1);
        c.showGame();
        Map<String, Double> shapleyValues = c.computeShapleyValues();
        boolean isShapleyInCore = c.isShapleyInCore(new ArrayList<>(shapleyValues.values()), new ArrayList<>(c.idMap.values()), c.nPlayers);
        System.out.println("Is Shapley value in the core? " + isShapleyInCore);
        boolean isCoreEmpty = isShapleyInCore ? false : isCoreEmpty(v1);
        System.out.println("Is core empty? " + isCoreEmpty);
    }

    public static boolean isCoreEmpty(double[] v1) {
        LinearProgramSolver solver = SolverFactory.newDefault();
        int nPlayers = (int) (Math.log(v1.length) / Math.log(2));
        double[] c = new double[nPlayers]; // Objective function
        Arrays.fill(c, 0.0);

        // Lower bound
        double[] lb = new double[nPlayers];
        Arrays.fill(lb, 0.0);

        // Create a new LinearProgram instance
        LinearProgram lp = new LinearProgram(c);
        lp.setMinProblem(true);

        // Add constraints
        for (int i = 0; i < Math.pow(2, nPlayers)-1; i++) {
            double[] A = new double[nPlayers];
            for (int j = 0; j < nPlayers; j++) {
                if ((i & (1 << j)) != 0) {
                    A[j] = 1.0;
                } else {
                    A[j] = 0.0;
                }
            }
            lp.addConstraint(new LinearBiggerThanEqualsConstraint(A, v1[i], "c" + i));
        }
        double[] A = new double[nPlayers];
        for(int i=0; i<nPlayers; i++)
            A[i] = 1.0;
        lp.addConstraint(new LinearEqualsConstraint(A, v1[(int) Math.pow(2, nPlayers)-1], "c" + (int) (Math.pow(2, nPlayers)-1)));
        lp.setLowerbound(lb);

        // Solve the linear program
        return solver.solve(lp) == null ? true : false;
    }

    private static double[] readVectorFromFile(String filename) {
        List<Double> list = new ArrayList<>();
        try (Scanner scanner = new Scanner(new File(filename))) {
            while (scanner.hasNextLine()) {
                list.add(Double.parseDouble(scanner.nextLine()));
            }
        } catch (FileNotFoundException e) {
        }
        double[] array = new double[list.size()];
        for (int i = 0; i < list.size(); i++) {
            array[i] = list.get(i);
        }
        return array;
    }

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

    public boolean isShapleyInCore(List<Double> shapleyValues, List<Integer> ids, int nPlayers) {
        List<Integer> players = new ArrayList<>();
        for (int i = 0; i < nPlayers; i++) {
            players.add(i);
        }
        Collections.reverse(ids);
        return checkAllSubsets(players, shapleyValues, ids, new ArrayList<>(), 0);
    }

    private boolean checkAllSubsets(List<Integer> players, List<Double> shapleyValues, List<Integer> ids, List<Integer> current, int index) {
        // If we have considered all players
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
        return v[index];
    }

    public Map<String, Double> computeShapleyValues() {
        System.out.println("*********** Computing Shapley Values ***********");
        Map<String, Double> shapleyValues = new TreeMap<>(); // TreeMap will sort the keys in natural order
        for (int i = 0; i < nPlayers; i++) {
            double shapleyValue = 0;
            for (int j = 0; j < Math.pow(2, nPlayers); j++) {
                if ((j & (1 << i)) != 0) {
                    // Compute the coalition without player i
                    int coalitionWithoutI = j - (1 << i);
                    // Compute the gain of player i when joining the coalition
                    double gain = v[j] - v[coalitionWithoutI];
                    // Compute the weight of the coalition
                    double weight = 1.0 / (nPlayers * binomialCoefficient(nPlayers - 1, Integer.bitCount(j) - 1));
                    // Update the Shapley value
                    shapleyValue += gain * weight;
                    System.out.println(returnShowSetString(coalitionWithoutI) + " (" + v[coalitionWithoutI] + ") -> "
                            + returnShowSetString(j)
                            + " (" + v[j] + ") ***** gain = " + gain + " (" + weight + ")");
                }
            }
            shapleyValue = Math.round(shapleyValue * 100.0) / 100.0;
            shapleyValues.put(ids[i], shapleyValue);
        }

        // Print Shapley values in alphabetical order based on the ids
        for (Map.Entry<String, Double> entry : shapleyValues.entrySet()) {
            System.out.println("Shapley Value for " + entry.getKey() + " : " + entry.getValue());
        }

        return shapleyValues;
    }

    public void setPlayersID() {
        int c = 64;
        ids = new String[nPlayers];
        for (int i = nPlayers - 1; i >= 0; i--) {
            c++;
            ids[i] = (String.valueOf((char) c));
        }
    }

    public void showGame() {
        System.out.println("*********** Coalitional Game ***********");
        for (int i = 0; i < v.length; i++) {
            showSet(i);
            System.out.println(" (" + v[i] + ")");
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
        boolean showPlayerID = true;
        // boolean showPlayerID = false;
        int power;
        System.out.print("{");
        int cnt = 0;
        for (int i = 0; i < nPlayers; i++) {
            power = nPlayers - (i + 1);
            if (showPlayerID) {
                if (inSet(i, v)) {
                    if (cnt > 0)
                        System.out.print(",");
                    cnt++;
                    System.out.print(ids[power]);
                }
            } else {
                if (cnt > 0)
                    System.out.print(",");
                cnt++;
                if (inSet(i, v))
                    System.out.print(1);
                else
                    System.out.print(0);
            }
        }
        System.out.print("}");
    }

    public boolean inSet(int i, long v) {
        int power;
        long vi;
        long div;
        long mod;
        power = nPlayers - (i + 1);
        vi = (long) Math.pow(2, power);
        div = v / vi;
        mod = div % 2;
        return (mod == 1);
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

    public void permutation(int j, int k, int iZero, long v0) {
        long value = 0;
        if (k == 0) {
            showSet(v0);
        } else {
            int op = 0;
            if (iZero < j)
                op = nPlayers - j;
            else
                op = nPlayers - j - 1;
            if (op == k) {
                for (int i = j; i < nPlayers; i++) {
                    if (i != iZero)
                        value += (long) Math.pow(2, nPlayers - (i + 1));
                }
                v0 = v0 + value;
                showSet(v0);
            } else {
                if (j != iZero)
                    permutation(j + 1, k - 1, iZero, v0 + (long) Math.pow(2, nPlayers - (j + 1)));
                permutation(j + 1, k, iZero, v0);
            }
        }
    }

    public static void showLP(LinearProgram lp) {
        System.out.println("*********** LINEAR PROGRAMMING PROBLEM ***********");
        String fs;
        if (lp.isMinProblem()) System.out.print("  minimize: ");
        else System.out.print("  maximize: ");
        double[] cf = lp.getC();
        for (int i = 0; i < cf.length; i++)
            if (cf[i] != 0) {
                fs = String.format(Locale.US, "%+7.1f", cf[i]);
                System.out.print(fs + "*x[" + i + "]");
            }
        System.out.println();
        System.out.print("subject to: ");
        ArrayList<Constraint> lcstr = lp.getConstraints();
        double aij;
        double[] ci = null;
        String str = null;
        for (int i = 0; i < lcstr.size(); i++) {
            if (lcstr.get(i) instanceof LinearSmallerThanEqualsConstraint) {
                str = " <= ";
                ci = ((LinearSmallerThanEqualsConstraint) lcstr.get(i)).getC();
            }
            if (lcstr.get(i) instanceof LinearBiggerThanEqualsConstraint) {
                str = " >= ";
                ci = ((LinearBiggerThanEqualsConstraint) lcstr.get(i)).getC();
            }
            if (lcstr.get(i) instanceof LinearEqualsConstraint) {
                str = " == ";
                ci = ((LinearEqualsConstraint) lcstr.get(i)).getC();
            }
            str = str + String.format(Locale.US, "%6.1f", lcstr.get(i).getRHS());
            if (i != 0) System.out.print("            ");
            for (int j = 0; j < lp.getDimension(); j++) {
                assert ci != null;
                aij = ci[j];
                if (aij != 0) {
                    fs = String.format(Locale.US, "%+7.1f", aij);
                    System.out.print(fs + "*x[" + j + "]");
                } else System.out.print("            ");
            }
            System.out.println(str);
        }
    }

}