// See https://aka.ms/new-console-template for more information
using System;
using System.Collections;
using System.Data;
using System.Reflection;
using System.Reflection.PortableExecutable;
using System.Text;
using System.Xml.Linq;
using Accord.Audio;
using Accord.Audio.Filters;
using Accord.Math;
using MathNet.Filtering.FIR;
using MathNet.Filtering.Windowing;
using MathNet.Numerics.LinearAlgebra;

Console.WriteLine("Hello, World!");

const double sampleRate = 400000;  //400KHZ
List<double> cutoffFrequency = new List<double>() { 5000, 10000, 15000, 20000, 25000 };  //50Hz

string method = "RCFilter";
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
    case "ZCI":
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
            if (!string.IsNullOrEmpty(data))
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
static double[] SubtractMean(double[] data1)
{
    double[] data = new double[data1.Length];
    double mean = data1.Average();
    for (int i = 0; i < data1.Length; i++)
    {
        data[i] = data1[i] - mean;
    }
    return data;
}

static int[] ZeroCrossings(double[] data)
{
    return Enumerable.Range(1, data.Length - 1)
        .Where(i => (data[i - 1] * data[i]) < 0)
        .ToArray();
}

static void SubtractAverages(ref double[] data, int start, int end, double average)
{
    for (int i = start; i < end; i++)
    {
        data[i] -= average;
    }
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
    Dictionary<Desc, double[]> keyValuePairs = new Dictionary<Desc, double[]>();
    int filterSize = 32;
    List<double> inputsig = GetInputSignal();
    double[] outputSignal1 = ApplyMedianFilter(inputsig.ToArray(), filterSize);

    inputsig = outputSignal1.ToList();

    keyValuePairs.Add(Desc.DecodedData, inputsig.ToArray());
    keyValuePairs.Add(Desc.MedianWithDecode, outputSignal1);

    double[] MeadianToMedian = ApplyMedianFilter(outputSignal1, filterSize);
    keyValuePairs.Add(Desc.MeadianToMeadian, MeadianToMedian);

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

    double[] inputSignal = lowpassOutput[0];
    keyValuePairs.Add(Desc.Meadian_LPFData, inputSignal);

    //double[] outputSignal = ApplyMedianFilter(inputSignal, filterSize);




    BitShiftFilter(lowpassOutput[0], keyValuePairs);




    WriteToCsv(keyValuePairs);
    wirteOutPutSignal(inputsig, lowpassOutput);

    static Dictionary<Desc, double[]> BitShiftFilter(double[] lowpassOutput, Dictionary<Desc, double[]> keyValuePairs)
    {

        double[] dataOut = lowpassOutput;

        // Subtract mean from dataOut
        double[] DatOut_Mean = SubtractMean(dataOut);

        // Find zero crossings
        int[] crossings = ZeroCrossings(DatOut_Mean);


        // Calculate bit-wise averages
        double[] averages = new double[crossings.Length];
        for (int idx = 0; idx < crossings.Length - 2; idx++)
        {
            double bitAMean = DatOut_Mean.Skip(crossings[idx]).Take(crossings[idx + 1] - crossings[idx]).Average();
            double bitBMean = DatOut_Mean.Skip(crossings[idx + 1]).Take(crossings[idx + 2] - crossings[idx + 1]).Average();
            averages[idx] = (bitAMean + bitBMean) / 2;
        }


        // Apply average to data
        double[] dataProcess = new double[dataOut.Length];
        Array.Copy(lowpassOutput, dataProcess, lowpassOutput.Length);

        for (int idx = 0; idx < crossings.Length - 1; idx++)
        {
            //if (idx == 0)
            //{
            //    SubtractAverages(ref dataProcess, 0, crossings[0], averages[0]);
            //}
            //else
            {
                SubtractAverages(ref dataProcess, crossings[idx], crossings[idx + 1], averages[idx]);
            }
        }
        //keyValuePairs.Add(Desc.MedianFilterData, dataOut);
        keyValuePairs.Add(Desc.OutputData, dataProcess);
        keyValuePairs.Add(Desc.DataOUT, DatOut_Mean);
        keyValuePairs.Add(Desc.Zci, crossings.Select(x => (double)x).ToArray());
        keyValuePairs.Add(Desc.Averages, averages);
        return keyValuePairs;
    }
}
static double[] ApplyMedianFilter(double[] inputSignal, int filterSize)
{
    int signalLength = inputSignal.Length;
    double[] outputSignal = new double[signalLength];

    for (int i = 0; i < signalLength; i++)
    {
        double[] neighbors = GetNeighbors(inputSignal, i, filterSize);
        outputSignal[i] = CalculateMedian(neighbors);
    }

    return outputSignal;
}

static double[] GetNeighbors(double[] signal, int index, int filterSize)
{
    int halfFilterSize = filterSize / 2;
    int startIndex = Math.Max(0, index - halfFilterSize);
    int endIndex = Math.Min(signal.Length - 1, index + halfFilterSize);

    double[] neighbors = new double[endIndex - startIndex + 1];
    Array.Copy(signal, startIndex, neighbors, 0, neighbors.Length);

    return neighbors;
}

static double CalculateMedian(double[] values)
{
    Array.Sort(values);

    int middle = values.Length / 2;

    if (values.Length % 2 == 0)
    {
        // If even number of elements, take average of the middle two
        return (values[middle - 1] + values[middle]) / 2;
    }
    else
    {
        // If odd number of elements, take the middle element
        return values[middle];
    }
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
void WritetOCSV(string Desc, double[] data)
{
    string outputFilePath = $"../../../../../BitShiftMech.csv";
    using (StreamWriter r = new StreamWriter(outputFilePath, true))
    {
        r.WriteLine();
        r.WriteLine(Desc);
        for (int i = 0; i < data.Length; i++)
        {
            r.WriteLine(data[i]);
        }
    }
}
static void WriteToCsv(Dictionary<Desc, double[]> dictionary)
{
    string filePath = $"../../../../../BitShiftMech1.csv";
    using (StreamWriter writer = new StreamWriter(filePath))
    {
        string Header = string.Join(",", dictionary.Keys.Select(key => key.ToString()));
        writer.WriteLine(Header);
        double maxCount = dictionary.Values.Max(values => values.Length);
        for (int i = 0; i < maxCount; i++)
        {
            string Body = string.Join(",", dictionary.Values.Select(values => values.Length > i ? values[i].ToString() : ""));
            writer.WriteLine($"{Body}");
        }
    }
}
static double FindMaxCount(Dictionary<Desc, double[]> dictionary)
{
    double maxCount = double.MinValue;

    foreach (var entry in dictionary)
    {
        double[] values = entry.Value;

        maxCount = Math.Max(maxCount, values.Length);
    }

    return maxCount;
}
public enum Desc
{
    DecodedData, MedianWithDecode, MeadianToMeadian, Meadian_LPFData, MedianFilterData, DataOUT, Zci, Averages, OutputData
};


