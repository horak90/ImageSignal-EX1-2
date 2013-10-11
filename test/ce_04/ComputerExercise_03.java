/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ce_04;

import java.io.ByteArrayInputStream;
import java.io.InputStream;
import lab.Signal;
import lab.ComplexSignal;
import utils.KeyboardIO;
import utils.Sound;

/**
 *
 * @author cfouard
 */
public class ComputerExercise_03 {

  /**
   * @param args the command line arguments
   */
  public static void main(String[] args) {
    int choice = -1;

    String menuMsg = "\n------------------------------------------------------------\n";
    menuMsg += " Main Menu:                                                 \n";
    menuMsg += "    1: Generate input signals                               \n";
    menuMsg += "    2: Test Discrete Fourier Transform                      \n";
    menuMsg += "    3: Test Inverse Discrete Fourier Transform              \n";
    menuMsg += "                                                            \n";
    menuMsg += "    4: Read and display a real    signal                    \n";
    menuMsg += "    5: Read and display a complex signal                    \n";
    menuMsg += "                                                            \n";
    menuMsg += "    9: Exit                                                 \n";
    menuMsg += " Please choose a number between 1 and 9                     \n";
    menuMsg += "------------------------------------------------------------\n\n";

    while (choice != 9) {
      choice = KeyboardIO.readInteger(menuMsg);

      String filename = null;
      String outputFilename = null;

      Signal signal = null;
      ComplexSignal complexSignal = null;
      String msg = null;

      switch (choice) {
        case 1:
          msg = "\n------------------------------------------------------------\n";
          msg += " Sub  Menu:                                                 \n";
          msg += "    1: Generate discrete delta signal                       \n";
          msg += "    2: Generate discrete square signal                      \n";
          msg += "    3: Generate complex sinc signal                         \n";
          msg += "    4: Generate sin                                         \n";
          msg += "    5: Generate cos                                         \n";
          msg += "    6: Generate a sound signal                              \n";
          msg += "                                                            \n";
          msg += "    9: Back to main menu                                    \n";
          msg += " Please choose a number between 1 and 9                     \n";
          msg += "------------------------------------------------------------\n\n";

          int choice2 = 0;
          while (choice2 != 9) {
            int nbSamples;
            int nonZeroSample;
            int firstNonZeroSample;
            int lastNonZeroSample;
            double a;
            int k;
            choice2 = KeyboardIO.readInteger(msg);
            switch (choice2) {
              case 1:
                nbSamples = KeyboardIO.readInteger("How many samples do you want to generate ?\n");
                nonZeroSample = KeyboardIO.readInteger("What is the index of the non zero sample?\n");

                signal = Signal.generateDelta(nonZeroSample, nbSamples);
                signal.display();

                // Saving signals
                if (KeyboardIO.readYesNo("Save this signal (please answer yes or no)?\n")) {
                  outputFilename = KeyboardIO.readString("Please enter the output file name to save delta\n");
                  signal.save(outputFilename);
                }
                break;

              case 2:
                nbSamples = KeyboardIO.readInteger("How many samples do you want to generate ?\n");
                firstNonZeroSample = KeyboardIO.readInteger("Index of the first non-zero sample: \n");
                lastNonZeroSample = KeyboardIO.readInteger("Index of the last  non-zero sample: \n");

                signal = Signal.generateRectangle(firstNonZeroSample, lastNonZeroSample, nbSamples);
                signal.display();

                // Saving signals
                if (KeyboardIO.readYesNo("Save this signal (please answer yes or no)?\n")) {
                  outputFilename = KeyboardIO.readString("Please enter the output file name to save rectangle signal\n");
                  signal.save(outputFilename);
                }
                break;

              case 3:
                nbSamples = KeyboardIO.readInteger("How many samples do you want to generate ?\n");
                int m = KeyboardIO.readInteger("Please enter the half width of the rectangle\n");
                complexSignal = Signal.generateComplexSinc(nbSamples, m);
                complexSignal.display();
                if (KeyboardIO.readYesNo("Do you want to display the real part only?\n")) {
                  complexSignal.displayReal();
                }
                if (KeyboardIO.readYesNo("Do you want to display the imaginary part only?\n")) {
                  complexSignal.displayImaginary();
                }
                if (KeyboardIO.readYesNo("Do you want to display the phase only?\n")) {
                  complexSignal.displayPhase();
                }
                if (KeyboardIO.readYesNo("Do you want to display the magnitude only?\n")) {
                  complexSignal.displayMagnitude();
                }

                if (KeyboardIO.readYesNo("Do you want to save the output signal?\n")) {
                  outputFilename = KeyboardIO.readString("Please enter the output file name.\n");
                  complexSignal.save(outputFilename);
                }
                break;

              case 4:
                nbSamples = KeyboardIO.readInteger("How many samples do you want to generate ?\n");
                a = KeyboardIO.readDouble("Please enter the amplitude");
                k = KeyboardIO.readInteger("Please enter frequency");
                signal = Signal.generateSin(nbSamples, a, k);
                signal.display();
                if (KeyboardIO.readYesNo("Do you want to save the output signal?\n")) {
                  outputFilename = KeyboardIO.readString("Please enter the output file name.\n");
                  signal.save(outputFilename);
                }
                break;

              case 5:
                nbSamples = KeyboardIO.readInteger("How many samples do you want to generate ?\n");
                a = KeyboardIO.readDouble("Please enter the amplitude");
                k = KeyboardIO.readInteger("Please enter frequency");
                signal = Signal.generateCos(nbSamples, a, k);
                signal.display();
                if (KeyboardIO.readYesNo("Do you want to save the output signal?\n")) {
                  outputFilename = KeyboardIO.readString("Please enter the output file name.\n");
                  signal.save(outputFilename);
                }
                break;
              case 6:
                filename = KeyboardIO.readString("Please enter the wav sound file name\n");
                Sound s = new Sound(filename);
                byte[] samples = s.getSamples();
                InputStream stream = new ByteArrayInputStream(samples);
                s.play(stream);

                nbSamples = KeyboardIO.readInteger("Please enter the number of sample you want (-1 for the whole signal, may be far too long for this tiny Java programm....)\n");
                signal = s.toSignal(nbSamples);
                signal.display();
                if (KeyboardIO.readYesNo("Do you want to save the output signal?\n")) {
                  outputFilename = KeyboardIO.readString("Please enter the output file name.\n");
                  signal.save(outputFilename);
                }
                break;


              case 9:
                break;
              default:
                System.out.println("I do not understand your choice, sorry...\n");
                break;
            }
          }
          break;


        case 2:
          filename = KeyboardIO.readString("Please enter the name of the signal file.\n");
          signal = new Signal(filename);
          complexSignal = signal.dft();
          complexSignal.display();
          if (KeyboardIO.readYesNo("Do you want to display the real part only?\n")) {
            complexSignal.displayReal();
          }
          if (KeyboardIO.readYesNo("Do you want to display the imaginary part only?\n")) {
            complexSignal.displayImaginary();
          }
          if (KeyboardIO.readYesNo("Do you want to display the phase only?\n")) {
            complexSignal.displayPhase();
          }
          if (KeyboardIO.readYesNo("Do you want to display the magnitude only?\n")) {
            complexSignal.displayMagnitude();
          }

          if (KeyboardIO.readYesNo("Do you want to save the output signal?\n")) {
            outputFilename = KeyboardIO.readString("Please enter the output file name.\n");
            complexSignal.save(outputFilename);
          }

          if (KeyboardIO.readYesNo("Do you want to save the magnitude?\n")) {
            outputFilename = KeyboardIO.readString("Please enter the output file name.\n");
            (complexSignal.getMagnitudeSignal()).save(outputFilename);
          }
          break;

        case 3:
          filename = KeyboardIO.readString("Where to find the complex signal?\n");
          complexSignal = new ComplexSignal(filename);

          signal = Signal.idft(complexSignal);

          signal.display();

          if (KeyboardIO.readYesNo("Do you want to save this signal?\n")) {
            outputFilename = KeyboardIO.readString("Please enter the output file name.\n");
            signal.save(outputFilename);
          }

          break;

        case 4:
          filename = KeyboardIO.readString("Please enter the name of the real signal file.\n");
          signal = new Signal(filename);
          signal.display();

          break;

        case 5:
          filename = KeyboardIO.readString("Please enter the name of the complex signal file.\n");
          complexSignal = new ComplexSignal(filename);
          complexSignal.display();
          if (KeyboardIO.readYesNo("Do you want to display the real part only?\n")) {
            complexSignal.displayReal();
          }
          if (KeyboardIO.readYesNo("Do you want to display the imaginary part only?\n")) {
            complexSignal.displayImaginary();
          }
          if (KeyboardIO.readYesNo("Do you want to display the phase only?\n")) {
            complexSignal.displayPhase();
          }
          if (KeyboardIO.readYesNo("Do you want to display the magnitude only?\n")) {
            complexSignal.displayMagnitude();
          }

          break;

        case 9:
          System.out.println("Good Bye.");
          break;
        default:
          System.out.println("The number you entered is not in the menu, sorry....");
          break;
      }
    }
  }
}
