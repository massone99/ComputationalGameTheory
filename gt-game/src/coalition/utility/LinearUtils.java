package coalition.utility;

import scpsolver.constraints.LinearEqualsConstraint;
import scpsolver.constraints.LinearSmallerThanEqualsConstraint;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;

public class LinearUtils {
    public LinearProgram setLP1(double[] minOrMaxFun, double[] constantTerms, double[][] coeffiecientsOfVariables, double[] lowerBound) {
        LinearProgram lp = new LinearProgram(minOrMaxFun);
        lp.setMinProblem(true);
        for (int i = 0; i < constantTerms.length; i++) {
            if (i == constantTerms.length - 1) {
                // Add the last constraint as an equality constraint
                lp.addConstraint(new LinearEqualsConstraint(coeffiecientsOfVariables[i], constantTerms[i], "coefficientsObjectiveFunction" + i));
            } else {
                // Add the other constraints as smaller than or equal to constraints
                lp.addConstraint(new LinearSmallerThanEqualsConstraint(coeffiecientsOfVariables[i], constantTerms[i], "coefficientsObjectiveFunction" + i));
            }
        }
        lp.setLowerbound(lowerBound);
        return lp;
    }

    public double[] solveLP(LinearProgram lp) {
        LinearProgramSolver solver = SolverFactory.newDefault();
        double[] x = solver.solve(lp);
        if (x == null) {
            throw new RuntimeException("No solution found");
        } else {
            return x;
        }
    }

    public void showSolution(LinearProgram lp, double[] x) {
        if (x == null) System.out.println("*********** NO SOLUTION FOUND ***********");
        else {
            System.out.println("*********** SOLUTION ***********");
            for (int i = 0; i < x.length; i++) System.out.println("x[" + i + "] = " + x[i]);
            System.out.println("f(x) = " + lp.evaluate(x));
        }
    }

    public double[][] utilityMatrixToEquations(double[][] utilityMatrix) {
        int numRows = utilityMatrix.length;
        int numCols = utilityMatrix[0].length;

        // Create a 2D array to store the matrix of equations for each column
        double[][] equationsMatrix = new double[numRows][numCols];

        // Iterate over each column in the utility matrix
        for (int j = 0; j < numCols; j++) {
            // Iterate over each row in the utility matrix
            for (int i = 0; i < numRows; i++) {
                // Store the value from the utility matrix in the equations matrix
                equationsMatrix[i][j] = utilityMatrix[i][j];
            }
        }

        return equationsMatrix;
    }

    private double[] getConstantTerms(int numConstantTerms) {
        double[] constantTerms = new double[numConstantTerms];
        for (int i = 0; i < numConstantTerms - 1; i++) {
            constantTerms[i] = 0;
        }
        constantTerms[numConstantTerms - 1] = 1;
        return constantTerms;
    }

    private double[][] addMinusOneColumnToMatrix(double[][] equationsMatrixP2) {
        int numRows = equationsMatrixP2.length;
        int numCols = equationsMatrixP2[0].length;

        // Create a new matrix with an additional column
        double[][] newMatrix = new double[numRows][numCols + 1];

        // Copy the values from the old matrix to the new one
        for (int i = 0; i < numRows; i++) {
            System.arraycopy(equationsMatrixP2[i], 0, newMatrix[i], 0, numCols);

            // Set the last element of each row to -1
            newMatrix[i][numCols] = -1;
        }

        return newMatrix;
    }

    public double[] getMinimizeFunctionCoefficients(double[][] equationMatrix) {
        int length = equationMatrix[0].length;
        double[] result = new double[length];
        for (int i = 0; i < length - 1; i++) {
            result[i] = 0;
        }
        result[length - 1] = 1;
        return result;
    }

    private void printEquationsMatrix(double[][] equationsMatrix) {
        System.out.println("Equations matrix:");
        for (double[] row : equationsMatrix) {
            for (double value : row) {
                System.out.print(value + " ");
            }
            System.out.println();
        }
    }
}
