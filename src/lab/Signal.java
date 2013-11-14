/*
 * Please modify this class and complete the following methods corresponding
 * to your computer exercise.
 */
package lab;

import java.text.DecimalFormat;
import java.util.ArrayList;

/**
 *
 * @author Horak Michael
 */
public class Signal extends GeneralSignal {

  public Signal() {
    super();
  }

  public Signal(String filename) {
    super(filename);
  }

  /* ************************************************************************* */
  /*                                                                           */
  /*         Computer Exercise number 01                                       */
  /*                                                                           */
  /* ************************************************************************* */
  /**
   * Computes and returns the mean value of the signal.
   * @return mean value of the signal.
   */
  public double getMean() 
  {
     int nbSamples = this.getNbSamples();
    double mean = 0.0;
    
    for (int i = 0; i < nbSamples; i++)
        mean += this.getValueOfIndex(i);

    mean /= nbSamples;

    return mean;
  }

  /**
   * Computes and returns the standard deviation (&#963;) of the signal.
   * @return Standard deviation (&#963;) of the signal.
   */
  public double getStandardDev() {
   int nbSamples = this.getNbSamples();
    double standardDev = 0.0;
    double mean = this.getMean();
    
    for (int i = 0; i < nbSamples; i++)
        standardDev += Math.pow(this.getValueOfIndex(i) - this.getMean(), 2);

    standardDev = Math.sqrt(standardDev / nbSamples);    
    return standardDev;
  }

  public double getVariance() {
    int nbSamples = this.getNbSamples();
    double variance = 0.0;
    // Write your code here

    return variance;
  }


  /**
   * Computes the histogramm of a signal for a given interval length.<br/>
   * This means that this function:
   *  - divides the value interval interval of the signal
   * ([getMin(), getMax()]) into N intervals of length {@link intervalLength}
   * N = ???.
   * <br/>
   * <ul>
   *    <li>for each interval of length intervalLength, counts the number of
   *        signal samples which value/ordinate is within the corresponding interval.
   *    </li>
   *    <li>stores this value at the ascissa corresponding to the center of the
   *        interval.
   *    </li>
   *
   * </ul>
   * @param intervalLength length of each interval
   * @return a signal containing the histogram of the current signal.
   */
  public Signal getHitogram(double intervalLength) {
      
    Signal histo = new Signal();
    if (intervalLength <= 0.0) {
      System.out.println("Error in computeHistogram method from Signal:\n The length of an interval must be strictly positive.\n");
      return null;
    }else
    {  
        DecimalFormat oneDForm = new DecimalFormat("#.#");        
        double minValue = Double.valueOf(oneDForm.format(this.getMin()));
        double maxValue = Double.valueOf(oneDForm.format(this.getMax()));    
        double intervalWidth = Double.valueOf(oneDForm.format((maxValue - minValue)/intervalLength));  
        double endOfInterval = Double.valueOf(oneDForm.format(minValue + (intervalLength * intervalWidth)-intervalWidth));
        
              
        for(double y=minValue; y <= endOfInterval; y += intervalWidth )
        {
            double tempInterval = Double.valueOf(oneDForm.format(y)); 
            int sumVal = 0;
            
            for(int x = 0; x < this.getNbSamples();x++)
            {         
                double sample = Double.valueOf(oneDForm.format(this.getValueOfIndex(x)));
                
                if(sample >= tempInterval && sample < tempInterval+intervalWidth)
                {                    
                    sumVal++;
                }else if(sample > endOfInterval && tempInterval==endOfInterval)
                {
                    sumVal++;
                }
            }             
            histo.addElement(tempInterval, sumVal);
        }           
         
     }

    return histo;
  }

   public static Signal createRandomSeries(int nbElements) {
    return createRandomSeries(0.0, 1.0, nbElements);
  }

  public static Signal createRandomSeries(double min, double max, int nbElements) {
    Signal resultSignal = new Signal();
    resultSignal.settName("Random Signal");

    double randomNumber;
    for (int i = 0; i < nbElements; i++) {
      randomNumber = min + Math.random() * (max - min);
      resultSignal.addElement(i, randomNumber);
    }
    return resultSignal;
  }

  public static Signal createGaussianNoiseSeries(int nbElements) {
    double val, random;
    Signal resultSignal = new Signal();
    resultSignal.settName("Gaussian Noise");
    for (int n = 0; n < nbElements; n++) {
      resultSignal.addElement(n, 0.0);
      for (int j = 0; j < 12; j++) {
        val = resultSignal.getValueOfIndex(n);
        random = Math.random();
        resultSignal.setValueOf(n, val + random);
      }
      val = resultSignal.getValueOfIndex(n);
      val -= 6;
      resultSignal.setValueOf(n, val);
    }
    return resultSignal;
  }

  

  /* ************************************************************************* */
  /*                                                                           */
  /*         Computer Exercise number 01   Team 01                             */
  /*                                                                           */
  /* ************************************************************************* */
  
  /**
   * Given the number of samples <i>nbSamples</i> (nbSamples = this.getNbSamples()) of the
   * current signal, decomposes the signal in <i>nbSamples</i> <b>impulse</b> component singals.
   * @return an ArrayList containing nbSamples impulse component signals.
   */
  public ArrayList<GeneralSignal> impulseDecomposition() {
    ArrayList<GeneralSignal> list = new ArrayList<GeneralSignal>();

    Signal signal = new Signal();
    list.add(signal);

    // Write your code here

    return list;
  }

  /**
   * Given the number of samples <i>nbSamples</i> (nbSamples = this.getNbSamples()) of the
   * current signal, decomposes the signal in <i>nbSamples</i> <b>step</b> component singals.
   * @return an ArrayList containing nbSamples step component signals.
   */
  public ArrayList<GeneralSignal> stepDecomposition() {
    ArrayList<GeneralSignal> list = new ArrayList<GeneralSignal>();
    // Write your code here

    return list;
  }

  /**
   * Decomposes the current signal into an <b>even signal</b> and an <b>odd signal</b>.
   * @return an ArrayList of 2 signals (one even component signal and one odd component signal).
   */
  public ArrayList<GeneralSignal> evenOddDecomposition() {
    ArrayList<GeneralSignal> list = new ArrayList<GeneralSignal>();
    // Write your code here

    return list;
  }

  /**
   * Decomposes the current signal into an <b>even samples signal</b> and an <b>odd samples signal</b>.
   * @return an ArrayList of 2 signals (one even component signal and one odd component signal).
   */
  public ArrayList<GeneralSignal> interlacedDecomposition() {
    ArrayList<GeneralSignal> list = new ArrayList<GeneralSignal>();
    // Write your code here

    return list;
  }

  /* ************************************************************************* */
  /*                                                                           */
  /*         Computer Exercise number 01   Team 02                             */
  /*                                                                           */


  
  /* ************************************************************************* */
  /*                                                                           */
  /*         Computer Exercise number 02                                       */
  /*                                                                           */
  /* ************************************************************************* */
  /**
   * Generates and returns the signal
   * x1[n] = sin(2 &Pi;n / 100)
   * @param nbSamples number of samples of the signal.
   * @return x1[n] = sin(2 &Pi;n / 100)
   */
  public static Signal generateX1(int nbSamples) {
    Signal x1 = new Signal();
    x1.settName("x1[n] = sin(2 Pi n / 100)");
    
    DecimalFormat oneDForm = new DecimalFormat("#.#");     
    
    for (int i = 0; i < nbSamples; i++)
    {
        double temp = Double.valueOf(oneDForm.format(Math.sin((2*(Math.PI)*i)/100)));        
        x1.addElement(i,temp);
    }   
    
    return x1;
  }

  /**
   * Generates and returns the signal
   * x2[n] =  4*exp(-(n-150)^2/300) - exp(-(n-150)^2/2500)
   * @param nbSamples number of samples of the signal.
   * @return x2[n] =  4*exp(-(n-150)^2/300) - exp(-(n-150)^2/2500)
   */
  public static Signal generateX2(int nbSamples) {
    Signal x2 = new Signal();
    x2.settName("x2[n] =  4*exp(-(n-150)^2/300) - exp(-(n-150)^2/2500)");    
    
    DecimalFormat oneDForm = new DecimalFormat("#.#");     
    
    for (int i = 0; i < nbSamples; i++)
    {
        
           double temp = Double.valueOf(oneDForm.format(
                   4 * (Math.exp((-1.0/300.0)*((i-150)^2)*1.0)) 
                   - (Math.exp((-1.0/2500.0)* ((i-150)^2)*1.0)))     );        
        x2.addElement(i,temp);        
    }   

    return x2;
  }

  /**
   * Generates and returns the signal
   * x3[n] =  1 for 240 &lt; n &gt; 300
   *         -2 for 299 &lt; n &gt; 380
   *          0  otherwise
   * @param nbSamples
   * @return x3
   */
  public static Signal generateX3(int nbSamples) {
    Signal x3 = new Signal();
    x3.settName("x3");    
    
    DecimalFormat oneDForm = new DecimalFormat("#.#");     
    
    for (int i = 0; i < nbSamples; i++)
    {
        double x;
        
        if(i > 240 && i < 300 )
        {
            x = 1;
        }else if(i > 299 && i < 380 )
        {
            x = -2;
        }else
        {
            x = 0;
        }
        
        double temp = Double.valueOf(oneDForm.format(x));        
        x3.addElement(i,temp);
    }   

    return x3;
  }

  /**
   *
   * Applies the following filter to current signal and returns the resulting signal:
   * y[0] = ??
   * For n from 1 to numberOfSamples-1
   * y[n] = 0.05*x[n] + 0.95*y[n-1]
   *
   */
  public Signal singlePolLowPassFilter() 
  {
    Signal result = new Signal();
    result.settName("Filtered signal");
    
    DecimalFormat oneDForm = new DecimalFormat("#.#");  
    
    int nbSamples = this.getNbSamples();
    result.addElement(0, 0.0);
    
    double y = 0.0;
    
    for(int i=1;i<nbSamples;i++)
    {        
        double temp = (0.05*this.getValueOfAbscissa(i) + 0.95 * y);
        y = temp;
                
        result.addElement(i,temp);  
        
    }
    
    
    return result;
  }

  /**
   * Linear convolve the current signal with the given kernel
   */
 public Signal linearConvolve(Signal kernel) {
    Signal result = new Signal();
    result.settName("Convolved signal");

    int n, k;
    int kMin, kMax;
    int kernelLength = kernel.getNbSamples();
    int signalLength = this.getNbSamples();
    int length = signalLength + kernelLength;
    double sum;
    double slope = this.getAbscissaOfIndex(1) - this.getAbscissaOfIndex(0);
    double C = this.getAbscissaOfIndex(0);

    for (n = 0; n < length - 1; n++) {
        sum = 0.0;
        kMin = (n >= kernelLength - 1) ? n - (kernelLength - 1) : 0;
        kMax = (n < signalLength - 1) ? n : signalLength - 1;

        for (k = kMin; k <= kMax; k++) 
        {
            sum += this.getValueOfIndex(k) * kernel.getValueOfIndex(n - k);
            System.out.print("n=" + n+ "kmin=" + kMin + "kmax=" + kMax + "k=" + k + "n-k" + (n-k) + "\n");
        }
       
        result.addElement(slope*n*1.0 + C, sum);
    }
   
    return result;
  }

  /**
   * Circular convolve the current signal with the given kernel
   */
/*
 
public Signal circularConvolve(Signal kernel)
{
    Signal result = new Signal();
    result.settName("Convolved signal");
   
    int circular = 0;
    double sum;
    
    for (int i = 0; i < this.getNbSamples(); i++) {
        sum = 0.0;
       
        for (int j = 0; j < kernel.getNbSamples(); j++) {
            circular = (i - j) < 0 ? (i - j + kernel.getNbSamples()) : i - j;
            sum += this.getValueOfIndex(j) * kernel.getValueOfIndex(circular);
        }
       
        result.addElement(this.getAbscissaOfIndex(i), sum);
    }   
    return result;

  }
*/
 
 
 
 
 
public Signal circularConvolve(Signal kernel) {
    Signal result = new Signal();
    result.settName("Convolved signal");
int end = kernel.getNbSamples()/2;
int start = (-1) * (kernel.getNbSamples()/2);
double sum;
int kshift;
for(int n = 0; n < this.getNbSamples(); n++){
sum = 0.0;
for (int k = start; k <= end; k++){
if((n-k >= 0)&& (n-k < getNbSamples()) ){
if (k < 0){
kshift = k + kernel.getNbSamples();
}else {
kshift = k;
}
//System.out.println("n=" + n + "kshift=" + kshift + "n-k" + (n-k) + "\n");
sum += this.getValueOfIndex(n-k)*kernel.getValueOfIndex(kshift);
}
}
result.addElement(n, sum);
}
return result;
}
  
  /* ************************************************************************* */
  /*                                                                           */
  /*         Computer Exercise number 02  Team 03                              */
  /*                                                                           */
  /* ************************************************************************* */
    // Write your additional functions here
  
  /* ************************************************************************* */
  /*                                                                           */
  /*         Computer Exercise number 03   Team 04                             */
  /*                                                                           */
  /* ************************************************************************* */
  // Write your additional functions here
  
  
  /* ************************************************************************* */
  /*                                                                           */
  /*         Computer Exercise number 03    (see Image.java)                   */
  /*                                                                           */
  /* ************************************************************************* */

  
  /* ************************************************************************* */
  /*                                                                           */
  /*         Computer Exercise number 03   Team 06   (see Image.java)          */
  /*                                                                           */
  /* ************************************************************************* */

  /* ************************************************************************* */
  /*                                                                           */
  /*         Computer Exercise number 03   Team 05                             */
  /*                                                                           */
  /* ************************************************************************* */
  public Signal stretchContrast(double newRangeMin, double newRangeMax) {
    Signal signal = new Signal();
    double a = newRangeMin;
    double b = newRangeMax;
    double c = this.getMin();
    double d = this.getMax();
    
    
    for(int i = 0; i < this.getNbSamples(); i++)
    {
        double s = (this.getValueOfIndex(i) - c) * ((b-a)/(d-c)) + a;
        signal.addElement(i, s);
    }   

    return signal;
  }


  /* ************************************************************************* */
  /*                                                                           */
  /*         Computer Exercise number 04                                       */
  /*                                                                           */
  /* ************************************************************************* */
  /**
   * Computes the Discrete Fourier Transform of the current signal
   * @return DFT of the current signal (as ComplexSignal)
   */
  public ComplexSignal dft() 
  {
    ComplexSignal result = new ComplexSignal();
    int nbSamples = this.getNbSamples();
    
    for(int n = 0; n < nbSamples; n++)
    {
        Complex cn= new Complex();
        
        for(int k = 0;k < nbSamples; k++)
        {
            double angle = -2*Math.PI*k*n/nbSamples;
            Complex ei = Complex.createFromPolar(1.0, angle);
            double xk = this.getValueOfIndex(k);
            ei.multiplyByReal(xk);
            cn.add(ei);
        }
        
        result.add(cn);
    }
return result;
}
  /**
   * Computes the Inverse Discrete Fourier Transform of the complex signal given
   * in parameters.<br/>
   * Note: this method should be in ComplexSignal, but for simplicity matters you
   * will not have to modify the class ComplexSignal, so you implement a static
   * method in Signal instead.
   * Note II: Compute the idft with complex number. At the end, take the real
   *   part of your result signal.
   * @param fourier input Fourier transform of the signal
   * @return the inverse Fourier transform of the input complex signal
   */
  public static Signal idft(ComplexSignal fourier) {
    ComplexSignal result = new ComplexSignal();
    int nbSamples = fourier.getNbSamples();
    // Write your code here

    return result.getRealSignal();
  }

  /**
   * Generates a &delta; signal
   * @param nonZeroSampleNumber number of the sample which value equals 1
   * @param numberOfSamples total number of samples of the result signal
   * @return Delta signal
   */
  public static Signal generateDelta(int nonZeroSampleNumber, int numberOfSamples) {
    Signal signal = new Signal();
    signal.settName("Delta_" + nonZeroSampleNumber);
    // Write your code here

    return signal;
  }

  /**
   * Generate a Rectangular signal of numberOfSamples samples:
   *  before firstNonZeroSample, and after lastNonZeroSample, sample values are 0.0
   *  between firstNonZeroSample and lastNonZeroSample included, sample values are 1.0
   * @param firstNonZeroSample index of the first non zero sample
   * @param lastNonZeroSample index of the last non zero sample
   * @param numberOfSamples number of samples of the signal
   */
  public static Signal generateRectangle(int firstNonZeroSample, int lastNonZeroSample, int numberOfSamples) {
    Signal signal = new Signal();
    signal.settName("Rectangle_" + firstNonZeroSample + "_" + lastNonZeroSample);
    // Write your code here

    return signal;
  }

  /**
   * Compute a complex signal where for each sample (of index n)
   * - the magnitude is
   *     - m if n=0
   *     - if n>0: sinc_m[n] = sin(2.pi.n/N . (m+1/2)) / sin(pi.n/N)
   *        (where N is the number of samples (numberOfSamples parameter))
   *
   * - the phase is 0
   *
   * @param numberOfSamples number of samples of the signal (N in the equation above)
   * @param m (as defined in the subject)
   * @return
   */
  public static ComplexSignal generateComplexSinc(int numberOfSamples, int m) {
    ComplexSignal signal = new ComplexSignal();
    // Write your code here

    return signal;
  }

  public static Signal generateXSignal(int numberOfSamples) {
    Signal signal = new Signal();
    signal.settName("(sin(2pi n 0.08)+2sin(2pi n 0.3))e^{-(n-200)^2 / 60^2}");
    // Write your code here

    return signal;
  }

  public static Signal generateSin(int nbSamples, double a, int k) {
    Signal signal = new Signal();
    // Write your code here

    return signal;
  }

  public static Signal generateCos(int nbSamples, double a, int k) {
    Signal signal = new Signal();
    // Write your code here

    return signal;
  }

}
