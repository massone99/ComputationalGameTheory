//package testing;
//
//import gametree.GameNode;
//import gametree.GameTree;
//
//import java.security.SecureRandom;
//import java.util.Iterator;
//
//import static lp.SolvingZeroSumGames.deserializeGameTree;
//
//public class TestingTree {
//
//    public static void main(String[] args) {
//        GameTree tree = deserializeGameTree();
//
//        // Init variables
//        GameNode rootNode = tree.getRootNode();
//
//        int n1 = rootNode.numberOfChildren();
//        int n2 = rootNode.getChildren().next().numberOfChildren();
//
//        String[] labelsP1 = new String[n1];
//        String[] labelsP2 = new String[n2];
//        int[][] U1 = new int[n1][n2];
//        int[][] U2 = new int[n1][n2];
//        GameNode childNode1, childNode2;
//        int j, i = 0;
//
//        System.out.println("ok");
//
//        // Find payoffs for each possible pure strategy profile
//        System.out.println("Assuming normal form game, see comment in the code.");
//        Iterator<GameNode> childrenNodes1 = rootNode.getChildren();
//
//        System.out.println("okok");
//
//        // Compute U1 and U2
//        while (childrenNodes1.hasNext()) {
//            childNode1 = childrenNodes1.next();
//            labelsP1[i] = childNode1.getLabel();
//            j = 0;
//            Iterator<GameNode> childrenNodes2 = childNode1.getChildren();
//            while (childrenNodes2.hasNext()) {
//                childNode2 = childrenNodes2.next();
//                if (i == 0)
//                    labelsP2[j] = childNode2.getLabel(); // Assuming normal form game (perfect information not
//                // supported)
//                U1[i][j] = childNode2.getPayoffP1();
//                U2[i][j] = childNode2.getPayoffP2();
//                j++;
//            }
//            i++;
//        }
//
//        System.out.println();
//    }
//}
