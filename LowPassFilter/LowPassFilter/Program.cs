// See https://aka.ms/new-console-template for more information
using System.Collections;
using System.Reflection;
using System.Reflection.PortableExecutable;
using System.Text;
using System.Xml.Linq;
using Accord.Audio;
using Accord.Audio.Filters;
using MathNet.Filtering.FIR;
using MathNet.Filtering.Windowing;

Console.WriteLine("Hello, World!");

const double sampleRate = 400000;  //400KHZ
List<double> cutoffFrequency = new List<double>() { 5000, 10000, 15000, 20000, 25000 };  //50Hz

string method = "Accord";
//string method = "MathNet";
//string method = "RCFilter";

switch (method)
{
    case "Accord":
        VerifyUsingAccord();
        break;
    case "MathNet":
        VerifyUsingMathNet();
        break;
    case "RCFilter":
        VerifyUsingRCFilter();
        break;
}

List<double> GetInputSignal()
{
    List<double> inputsig = new List<double>();
    //string[] allLines = File.ReadAllLines("../../../../../inputSignalRawData.csv");
    string[] allLines = File.ReadAllLines("../../../../../inputSignalDecodeData.csv");
    if (allLines != null)
    {
        for (long i = 0; i < allLines.Length; i++)
        {
            string data = allLines[i];
            inputsig.Add(double.Parse(data));

        }
    }
    return inputsig;
}


void VerifyUsingMathNet()
{
    List<double> inputsig = GetInputSignal();
    List<double[]> lowpassOutput = new List<double[]>();
    for (int i = 0; i < cutoffFrequency.Count; i++)
    {
        var lowPass = OnlineFirFilter.CreateLowpass(MathNet.Filtering.ImpulseResponse.Finite,
            sampleRate, cutoffFrequency[i]);
        double[] outputArray = lowPass.ProcessSamples(inputsig.ToArray());
        lowpassOutput.Add(outputArray);
    }
    wirteOutPutSignal(inputsig, lowpassOutput);
}

void VerifyUsingAccord()
{
    //// Apply low-pass filter using Accord.NET
    List<double> inputsig = GetInputSignal();
    List<double[]> lowpassOutput = new List<double[]>();
    for (int i = 0; i < cutoffFrequency.Count; i++)
    {
        Signal inputAccordSignal = Signal.FromArray(inputsig.ToArray(), (int)sampleRate);
        LowPassFilter lowPassFilter = new LowPassFilter(cutoffFrequency[i], sampleRate);
        Signal outputAccordSignal = lowPassFilter.Apply(inputAccordSignal);
        byte[] outputSignal = outputAccordSignal.RawData;
        double[] outputArray = outputSignal.Select(i => (double)i).ToArray();
        lowpassOutput.Add(outputArray);
    }
    wirteOutPutSignal(inputsig, lowpassOutput);
}

void VerifyUsingRCFilter()
{
    //function lowpass(real[1..n] x, real dt, real RC)
    //    var real[1..n] y
    //    var real α := dt / (RC + dt)
    //    y[1] := α * x[1]
    //    for i from 2 to n
    //        y[i] := α * x[i] + (1 - α) * y[i - 1]
    //    return y
    List<double> inputsig = GetInputSignal();
    List<double> outSig = new List<double>();
    List<double[]> lowpassOutput = new List<double[]>();
    for (int j = 0; j < cutoffFrequency.Count; j++)
    {
        double RC = 1 / (2 * Math.PI * cutoffFrequency[0]);
        double dt = 1 / sampleRate;
        float alfa = (float)(dt / (dt + RC));

        outSig.Add(alfa * inputsig[0]);
        for (int i = 1; i < inputsig.Count; i++)
        {
            outSig.Add(alfa * inputsig[i] + (1 - alfa) * outSig[i - 1]);
        }
        lowpassOutput.Add(outSig.ToArray());
    }
    wirteOutPutSignal(inputsig, lowpassOutput);
}

void wirteOutPutSignal(List<double> inputsig, List<double[]> outSig)
{
    string outputFilePath = $"../../../../../RClowpassFilter_usingDecData_{DateTime.Now.ToString("MM-dd-yyyy_hh-mm-ss-tt")}.csv";
    using (StreamWriter r = new StreamWriter(outputFilePath, true))
    {
        r.WriteLine($"Input Signal,Filtered Signal {cutoffFrequency[0]}, Filtered Signal {cutoffFrequency[1]}," +
            $"Filtered Signal {cutoffFrequency[2]}, Filtered Signal {cutoffFrequency[3]}, " +
            $"Filtered Signal {cutoffFrequency[4]}");
        for (int i = 0; i < inputsig.Count; i++)
        {
            r.WriteLine($"{inputsig[i]:F3},{outSig[0][i]:F3}," +
                $"{outSig[1][i]:F3},{outSig[2][i]:F3},{outSig[3][i]:F3}, {outSig[4][i]:F3}");
        }
    }
}


