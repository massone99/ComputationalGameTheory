package lp;

import java.util.ArrayList;
import java.util.Locale;

import scpsolver.constraints.*;
import scpsolver.lpsolver.*;
import scpsolver.problems.*;


public class LinearProgrammingFuncs {

    public LinearProgrammingFuncs() {
    }


    public static LinearProgram setLP1() {
        double[] coefficientsObjectiveFunction = {150.0, 175.0};
        double[] constantTerms = {77.0, 80.0, 9.0, 6.0};
        double[][] coeffiecientsOfVariables = {
                {7.0, 11.0},
                {10.0, 8.0},
                {1.0, 0.0},
                {0.0, 1.0},
        };
        double[] lowerBound = {0.0, 0.0};
        LinearProgram lp = new LinearProgram(coefficientsObjectiveFunction);
        lp.setMinProblem(false);
        for (int i = 0; i < constantTerms.length; i++)
            lp.addConstraint(new LinearSmallerThanEqualsConstraint(coeffiecientsOfVariables[i],
                    constantTerms[i],
                    "coefficientsObjectiveFunction" + i));
        lp.setLowerbound(lowerBound);
        return lp;
    }


    public static double[] solveLP(LinearProgram lp) {
        LinearProgramSolver solver = SolverFactory.newDefault();
        double[] x = solver.solve(lp);
        if (x == null) {
            throw new RuntimeException("No solution found");
        } else {
            return x;
        }
    }

    public static void showSolution(LinearProgram lp, double[] x) {
        if (x == null) System.out.println("*********** NO SOLUTION FOUND ***********");
        else {
            System.out.println("*********** SOLUTION ***********");
            for (int i = 0; i < x.length; i++) System.out.println("x[" + i + "] = " + x[i]);
            System.out.println("f(x) = " + lp.evaluate(x));
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
        System.out.println("");
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
                aij = ci[j];
                if (aij != 0) {
                    fs = String.format(Locale.US, "%+7.1f", aij);
                    System.out.print(fs + "*x[" + j + "]");
                } else System.out.print("            ");
            }
            System.out.println(str);
        }
    }


    public static void main(String[] args) {
        LinearProgram lp = setLP1();
        showLP(lp);
        double[] x = solveLP(lp);
        showSolution(lp, x);
    }

}
