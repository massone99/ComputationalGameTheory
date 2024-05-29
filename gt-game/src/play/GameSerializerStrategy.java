package play;

import gametree.GameTree;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;


public class GameSerializerStrategy extends Strategy {

    public void execute() throws InterruptedException {

        while (!this.isTreeKnown()) {
            System.err.println("Waiting for game tree to become available.");
            Thread.sleep(1000);
        }

        serializeGameTree(tree, "TT.ser");

        while (true) {

            PlayStrategy myStrategy = this.getStrategyRequest();
            if (myStrategy == null) // Game was terminated by an outside event
                break;

            boolean playComplete = false; // Strategy was not valid

            while (!playComplete) {
            }

        }
    }

    public void serializeGameTree(GameTree tree, String filePath) {
        try {
            FileOutputStream fileOut = new FileOutputStream(filePath);
            ObjectOutputStream out = new ObjectOutputStream(fileOut);
            out.writeObject(tree);
            out.close();
            fileOut.close();
        } catch (IOException i) {
            i.printStackTrace();
        }
    }


}