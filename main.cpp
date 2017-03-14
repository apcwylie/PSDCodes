//Code written to process the data from a Scionix VS-1161-10 detector in various environments.
//The pulses recorded are here submitted to pulse shape analysis (PSA) to separate out the
//neutron like particles from the rest of the particles present in the detector output.
//All times are in seconds.
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>

#define M_PI 3.14159265358979323846
//The dimensions of the active component of the detector
#define EJ426DETX 50 //centimetres
#define EJ426DETY 5.1 //centimetres
#define EJ426DETZ 0.032 //centimetres
#define DETY 5.4 //cm
#define DETZ 2.35 //cm
#define DETAREAERROR 0.5 //cm^2
#define TIMEERR 0.5 //seconds


using namespace std;

//------------------------------------------------Some utilities first--------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------

//Method to return the value in a vector<double> furthest from 0.
double maxModVal(vector<double> input){
    if (input.size() == 0){
        cout<< "No elements" << endl;
        return 0;
    }
    double output = input[0];
    for (int i=0; i<input.size();++i){
        if (abs(input[i])>abs(output)){
            output = input[i];
        }
    }
    return output;
}

//Method to return the modulus of the value in a vector<double> furthest from 0.
double modMaxModVal(vector<double> input){
    if (input.size() == 0){
        cout<< "No elements" << endl;
        return 0;
    }
    double output = input[0];
    for (int i=0; i<input.size();++i){
        if (abs(input[i])>abs(output)){
            output = abs(input[i]);
        }
    }
    return output;
}

//Method to count the number of waves in a file.
void numWaves(string inFileName, int wSize){
    int counter, time, numWaves;
    fstream f_in;
    vector<double> wave;
    double height;
    f_in.open(inFileName.c_str(),std::fstream::in);
    if(!f_in){
        cout<< " not found in numWaves with filename: " + inFileName << endl;
    }
    numWaves = 0;
    counter = 0;
    //Start reading in values.
    f_in >> time >> height;
    while(f_in){
        if(counter<wSize){
            wave.push_back(height); // Same wave
            counter++;
        }
        else {
            numWaves++;
            wave.clear();
            counter = 0;
        }
        f_in >> time >> height;
    }
    f_in.close();
    cout<<"The number of waves in "<<inFileName<<" is: "<<numWaves<<endl;
}

//Method to print the figure of merit and its error from input peak separation
//and peak widths with errors.
void FoM(double X, double dX, double W_a, double dW_a, double W_b, double dW_b){
    double figure = X/(W_a+W_b);
    double dFoM = sqrt(dX*dX/((W_a+W_b)*(W_a+W_b))
                        +(X*dW_a/((W_a+W_b)*(W_a+W_b)))*(X*dW_a/((W_a+W_b)*(W_a+W_b)))
                          + (X*dW_b/((W_a+W_b)*(W_a+W_b)))*(X*dW_b/((W_a+W_b)*(W_a+W_b))));
    cout << "The figure of merit for the inputs is: " << figure<< " with an error of: " << dFoM <<endl;
}


//Method to print the first 10 waveforms in a file to a txt file.
void firstTen(string inFileName, string outFileName, int wSize){
    //Clear output file and set up variables.
    ofstream f_outClear;
    f_outClear.open(outFileName, std::ofstream::out | std::ofstream::trunc);
    f_outClear.close();
    int counter, numWaves, time;
    double height;
    counter = numWaves = 0;
    fstream f_in;
    vector<double> wave1;
    f_in.open(inFileName.c_str(),std::fstream::in);
    ofstream f_out(outFileName, ios::out | ios::app);
    if(!f_in){
        cout<< " not found in firstTen with filename: " + inFileName<< endl;
    }
    //Start reading in values.
    f_in >>time>> height;
    while(f_in){
        if(counter<wSize){
            wave1.push_back(height); // Same wave
            counter++;
        }
        else {
            numWaves++;
            //Save values
            if (numWaves<10) {
                if (f_out.is_open()) {
                    for (int i = 0; i < wave1.size(); ++i) {
                        f_out <<i <<" "<< wave1[i]<<endl;
                    }
                } else {
                    cout << "Unable to open file " << endl;
                }
            }else{
                break;
            }
            wave1.clear();
            counter = 0;
        }
        f_in>>time>>height;
    }
    f_in.close();
    f_out.close();
    cout<<"                       firstTen Completed                    "<<endl;

}

//Method to adjust the values in a file by their baseline to zero. The first baseLEnd values are
//averaged and then subtracted from the whole wave.
void baselineAdjust(string inFileName, string outFileName, int wSize, int baseLEnd){
    //Clear output file and set up variables.
    ofstream f_outClear;
    f_outClear.open(outFileName, std::ofstream::out | std::ofstream::trunc);
    f_outClear.close();
    int counter, time;
    double height, basel;
    counter = basel = 0;
    fstream f_in;
    vector<double> wave;
    f_in.open(inFileName.c_str(),std::fstream::in);
    ofstream f_out(outFileName, ios::out | ios::app);
    if(!f_in){
        cout<< " not found in baselineAdjust with filename: " + inFileName<< endl;
    }
    //Start reading in values.
    f_in >> time>> height;
    while(f_in){
        if(counter<wSize){
            wave.push_back(height); // Same wave
            counter++;
        }
        else {
            //calculate baseline.
            basel = 0;
            for (int i = 0; i < baseLEnd; ++i) {
                basel += wave[i]/baseLEnd;
            }
            //Subtract baseline
            for (int i = 0; i < wave.size(); ++i) {
                wave[i] -= basel;
            }
            //Save values
            if (f_out.is_open()) {
                for (int i = 0; i < wave.size(); ++i) {
                    f_out << i << " "<<wave[i] << endl;
                }
            } else {
                cout << "Unable to open file " << endl;
            }
            wave.clear();
            counter = 0;
        }
        f_in>>time>>height;
    }
    f_in.close();
    f_out.close();
    cout<<"                       baselineAdjust Completed                    "<<endl;

}

//-----------------------------------------------------PSA METHODS------------------------------------------------------
//Forms the majority of the code, several methods are attempted here, the most notable of which being thewidths method
//which shows promising results in discriminating the heavy particle signature (neutrons) from the other signals
// (thought to be gammas mostly)
//----------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------Tail vs Peak Integral-----------------------------------------------
//----------------------------------------------------------------------------------------------------------------------


//Giant method to do all of the peak/tail PSA all at once.
//PEAK IS FIRST COLUMN, TAIL IS SECOND IN OUTPUT FILE
void peakTailIntegrate(string inFileName, string outFileName, int wSize, int baseLEnd, int peakXValue, int tailEndXVal){
    //Clear output file and set up variables.
    ofstream f_outClear;
    f_outClear.open(outFileName, std::ofstream::out | std::ofstream::trunc);
    f_outClear.close();
    int  pulserNo, counter,peakXVal;
    double height;
    pulserNo = counter = 0;
    fstream f_in;
    vector<double> wave;
    double peak,tail,basel, maxVal,time;
    f_in.open(inFileName.c_str(),std::fstream::in);
    ofstream f_out(outFileName, ios::out | ios::app);
    if(!f_in){
        cout<< " not found in peakTailIntegrate with filename: " + inFileName<< endl;
    }
    //Start reading in values.
    f_in >> time>>height;
    while(f_in){
        if(counter<wSize){
            wave.push_back(height); // Same wave
            counter++;
        }
        else {
            //remove pulsers.
            peak = tail = basel = 0;
            for (int i=0;i<baseLEnd;++i){
                basel += wave[i]/baseLEnd;
            }
            for (int i=0;i<wSize;++i){
                wave[i] -= basel;
            }
            maxVal = maxModVal(wave);
            //Integrate.
            for (int i = 0; i < wave.size(); ++i) {
                if (i < peakXValue){
                    peak += wave[i];
                }
                else if ((i > peakXValue) && (i < tailEndXVal)){
                    tail += wave[i];
                }
            }
            //Save values
            if (f_out.is_open()) {
                f_out << peak << " " << tail << endl;
            } else {
                cout << "Unable to open file " << endl;
            }
            wave.clear();
            counter = 0;
        }
        f_in>>time>>height;
    }
    f_in.close();
    f_out.close();
    cout<<"                       peakTailIntegrate Completed                    "<<endl;

}

//-------------------------------------------Risetime vs Peak Height----------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------


//Method to calculate integral rise time vs amplitude for an input file
//PEAK VALUE IS FIRST COLUMN, INTEGRAL RISETIME IS SECOND COLUMN).
void IntegralRisetimeVsAmplitude(string inFileName, string outFileName, double lowThresh, double highThresh, int wSize, int baseLEnd){
    //Clear output file and set up variables.
    ofstream f_outClear;
    f_outClear.open(outFileName, std::ofstream::out | std::ofstream::trunc);
    f_outClear.close();
    int counter, lowTime, highTime, risetime;
    fstream f_in;
    vector<double> wave;
    double basel, totalIntegral, peak, accumulate, height;
    f_in.open(inFileName.c_str(),std::fstream::in);
    ofstream f_out(outFileName, ios::out | ios::app);
    if(!f_in){
        cout<< " not found in IntegralRisetimeVsAmplitude with filename: " + inFileName<< endl;
    }
    counter = 0;
    //Start reading in values.
    f_in >> height;
    while(f_in){
        if(counter<wSize){
            wave.push_back(height); // Same wave
            counter++;
        }
        else {
            //Calculate baseline.
            basel = totalIntegral = accumulate = 0;
            for (int i=0; i < baseLEnd; ++i) {
                basel += wave[i]/baseLEnd;
            }
            //Subtract baseline an integrate.
            for (int i=0; i<wSize; i++){
                wave[i] -= basel;
                totalIntegral += wave[i];
            }
            peak = maxModVal(wave);
            for (int i=0; i<wave.size();++i){
                accumulate += wave[i];
                if (accumulate > lowThresh*totalIntegral){
                    lowTime = i;
                    break;
                }
            }
            accumulate = 0;
            for (int i=0; i<wSize;++i){
                accumulate += wave[i];
                if (accumulate > highThresh*totalIntegral){
                    highTime = i;
                    break;
                }
            }
            risetime = highTime - lowTime;
            //Save values.
            if (f_out.is_open()) {
                f_out << peak << " " << risetime << endl;
            } else {
                cout << "Unable to open file " << endl;
            }
            wave.clear();
            counter = 0;
        }
        f_in>>height;
    }
    f_in.close();
    f_out.close();
    cout<<"                       IntegralRisetimeVsAmplitude Completed                    "<<endl;

}

//--------------------------------------------------------------Widths--------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------

//Method to calculate the width of the pulse for a given fraction of its height, such as the full width half maximum
//This one works for an input file that is a list of wave heights of size WSIZE
void Widths(string inFileName, string outFileName, double threshold, int wSize, int baseLEnd){
    //Clear output file and set up variables.
    ofstream f_outClear;
    f_outClear.open(outFileName, std::ofstream::out | std::ofstream::trunc);
    f_outClear.close();
    int counter, lowTime, highTime, width, time;
    fstream f_in;
    vector<double> wave;
    double maxVal, height, basel;
    f_in.open(inFileName.c_str(),std::fstream::in);
    ofstream f_out(outFileName, ios::out | ios::app);
    if(!f_in){
        cout<< " not found in Widths with filename: " + inFileName<< endl;
    }
    counter = 0;
    //Start reading in values.
    f_in >> time >> height;
    while(f_in){
        if(counter<wSize){
            wave.push_back(height); // Same wave
            counter++;
        }
        else {
            basel = 0.0;
            //Subtract the baseline (For LUNA results this looks like around 2244?).
            for(int i=0;i<baseLEnd;i++){
                basel+=wave[i];
            }
            basel = basel/baseLEnd;
            for(int i=0;i<wave.size();++i){
                wave[i]-=basel;
            }
            //Find maxVal for the wave.
            maxVal = maxModVal(wave);
            for (int i=0; i<wave.size();++i){
                if (abs(wave[i]) > threshold*abs(maxVal)){
                    lowTime = i;
                    break;
                }
            }
            //Find the width.
            for (int i=wSize; i>0;--i){
                if (abs(wave[i]) > threshold*abs(maxVal)){
                    highTime = i;
                    break;
                }
            }
            width = highTime - lowTime;
            //eliminate the noise cases with widths of 3999 or similar and output them.
            if ((width < 0.8*wSize)&&(width>0.0)){
                if (f_out.is_open()){
                    f_out << width << endl;
                } else {
                    cout << "Unable to open file: " + outFileName<< endl;
                }
            }
            wave.clear();
            counter = 0;
        }
        f_in >> time >> height;
    }
    f_in.close();
    f_out.close();
    cout<<"                       Widths Completed                    "<<endl;

}


//Method to bin the data (in bin sizes of 1) from one of the output files from the Widths method.
//This method normalises the data to
void widthBinTimeNormalised(string inFileName, string outFileName, double time, double binSize, int wSize){
    //Clear output file and set up variables.
    ofstream f_outClear;
    f_outClear.open(outFileName, std::ofstream::out | std::ofstream::trunc);
    f_outClear.close();
    double width;
    fstream f_in;
    vector<double> widthBinVals((int)(wSize/binSize));
    f_in.open(inFileName.c_str(),std::fstream::in);
    if(!f_in){
        cout<< " not found in widthBinTimeNormalised with filename: " + inFileName<< endl;
    }
    //read in values.
    f_in >> width;
    while(f_in){
        int index = round(width/binSize);
        widthBinVals[index]++;
        f_in >> width;
    }
    f_in.close();
    ofstream f_out(outFileName, ios::out | ios::app);

    for(int i=0;i<widthBinVals.size();++i){
        widthBinVals[i]/=time;
        if (f_out.is_open()) {
            f_out << i * binSize << " " << widthBinVals[i] << endl;
        }else{
            cout <<"Unable to open file: " + outFileName <<endl;
        }
    }
    f_out.close();
    cout<<"                       widthBinTimeNormalised Completed                    "<<endl;

}


//-------------------------------------------------Total Integral vs Width----------------------------------------------
//----------------------------------------------------------------------------------------------------------------------

//Method to calculate total integral and width of a set of waveforms.
//It also gets rid of any extremely wide waves of width greater than 80% wave size
//WIDTH IS FIRST COLUMN, TOTAL INTEGRAL IS SECOND IN OUTPUT FILE
void totalIntVsWidth(string inFileName, string outFileName, double threshold, int wSize, int baseLEnd, int wStart, int wEnd){
    //Clear output file and set up variables.
    ofstream f_outClear;
    f_outClear.open(outFileName, std::ofstream::out | std::ofstream::trunc);
    f_outClear.close();
    int counter, lowTime, highTime, width, time;
    fstream f_in;
    vector<double> wave;
    double maxVal, height, basel, totalInt;
    f_in.open(inFileName.c_str(),std::fstream::in);
    ofstream f_out(outFileName, ios::out | ios::app);
    if(!f_in){
        cout<< " not found in totalIntVsWidth with filename: " + inFileName<< endl;
    }
    counter = 0;
    //Start reading in values.
    f_in >> time >> height;
    while(f_in){
        if(counter<wSize){
            wave.push_back(height); // Same wave
            counter++;
        }
        else {
            basel = 0.0;
            //Subtract the baseline (For LUNA results this looks like around 2244?).
            for(int i=0;i<baseLEnd;i++){
                basel+=wave[i]/baseLEnd;
            }
            for(int i=0;i<wave.size();++i){
                wave[i]-=basel;
            }
            //Find maxVal for the wave.
            maxVal = maxModVal(wave);
            for (int i=0; i<wave.size();++i){
                if (abs(wave[i]) > threshold*abs(maxVal)){
                    lowTime = i;
                    break;
                }
            }
            //Find the width.
            for (int i=wSize; i>0;--i){
                if (abs(wave[i]) > threshold*abs(maxVal)){
                    highTime = i;
                    break;
                }
            }
            width = highTime - lowTime;
            //eliminate the noise cases with widths of 3999 or similar and output them.
            if ((width < 0.8*wSize)&&(width>0.0)){
                //Integral bit.
                totalInt = 0;
                for (int i=wStart; i<wEnd; ++i) {
                    totalInt += wave[i];
                }

                //Save values
                if (f_out.is_open()) {
                    f_out << width << " " << totalInt << endl;
                } else {
                    cout << "Unable to open file: " + outFileName << endl;
                }
            }
            wave.clear();
            counter = 0;
        }
        f_in >> time >> height;
    }
    f_in.close();
    f_out.close();
    cout<<"                       totalIntVsWidth Completed                    "<<endl;

}


//After the above baseline adjustment has been made to the file, this method will calculate the total int vs width.
void totalIntVsWidthPostBaselineAdjusted(string inFileName, string outFileName, double threshold, int wSize, int wStart, int wEnd){
    //Clear output file and set up variables.
    ofstream f_outClear;
    f_outClear.open(outFileName, std::ofstream::out | std::ofstream::trunc);
    f_outClear.close();

    int counter, lowTime, highTime, width, time;
    vector<double> wave;
    double maxVal, height, totalInt;

    fstream f_in;
    f_in.open(inFileName.c_str(),std::fstream::in);
    ofstream f_out(outFileName, ios::out | ios::app);
    if(!f_in){
        cout<< " not found in totalIntVsWidthPostBaselineAdjusted with filename: " + inFileName<< endl;
    }
    counter = 0;
    //Start reading in values.
    f_in >> time >> height;
    while(f_in){
        if(counter<wSize){
            wave.push_back(height); // Same wave
            counter++;
        }
        else {
            //Find maxVal for the wave.
            maxVal = maxModVal(wave);
            for (int i=0; i<wave.size();++i){
                if (abs(wave[i]) > threshold*abs(maxVal)){
                    lowTime = i;
                    break;
                }
            }
            //Find the width.
            for (int i=wSize; i>0;--i){
                if (abs(wave[i]) > threshold*abs(maxVal)){
                    highTime = i;
                    break;
                }
            }
            width = highTime - lowTime;
            //eliminate the noise cases with widths of 3999 or similar and output them.
            if ((width < 0.8*wSize)&&(width>0.0)){
                //Integral bit.
                totalInt = 0;
                for (int i=wStart; i<wEnd; ++i) {
                    totalInt += wave[i];
                }

                //Save values
                if (f_out.is_open()) {
                    f_out << width << " " << totalInt << endl;
                } else {
                    cout << "Unable to open file: " + outFileName << endl;
                }
            }
            wave.clear();
            counter = 0;
        }
        f_in >> time >> height;
    }
    f_in.close();
    f_out.close();
    cout<<"                       totalIntVsWidthPostBaselineAdjusted Completed                    "<<endl;

}

void totalIntPostBaselineAdjusted(string inFileName, string outFileName, int wSize, int wStart, int wEnd){
    //Clear output file and set up variables.
    ofstream f_outClear;
    f_outClear.open(outFileName, std::ofstream::out | std::ofstream::trunc);
    f_outClear.close();

    int counter, lowTime, highTime, time;
    vector<double> wave;
    double  height, totalInt;

    fstream f_in;
    f_in.open(inFileName.c_str(),std::fstream::in);
    ofstream f_out(outFileName, ios::out | ios::app);
    if(!f_in){
        cout<< " not found in totalIntPostBaselineAdjusted with filename: " + inFileName << endl;
    }
    counter = 0;
    //Start reading in values.
    f_in >> time >> height;
    while(f_in){
        if(counter<wSize){
            wave.push_back(height); // Same wave
            counter++;
        }
        else {
            //Find maxVal for the wave.
            //Integral bit.
            totalInt = 0;
            for (int i=wStart; i<wEnd; ++i) {
                totalInt += wave[i];
            }
            //Save values
            if (f_out.is_open()) {
                f_out << totalInt << endl;
            } else {
                cout << "Unable to open file: " + outFileName << endl;
            }
            wave.clear();
            counter = 0;
        }
        f_in >> time >> height;
    }
    f_in.close();
    f_out.close();
    cout<<"                       totalIntPostBaselineAdjusted Completed                    "<<endl;

}

//-----------------------------------------------Pulse Gradient Analysis------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
//Methods to perform the pulse gradient analysis detailed in Radiation Detection and Measurement (G.F.Knoll, 1989) comparing
//the (baseline adjusted) amplitude to a sample value.
void PGA(string inFileName, string outFileName, int sampleNo, int wSize, int baseLEnd ){
    ofstream f_outClear;
    f_outClear.open(outFileName, std::ofstream::out | std::ofstream::trunc);
    f_outClear.close();
    int counter, time;
    fstream f_in;
    vector<double> wave;
    double amplitudeVal, height, basel, sampleVal, PGAVal;
    f_in.open(inFileName.c_str(),std::fstream::in);
    ofstream f_out(outFileName, ios::out | ios::app);
    if(!f_in){
        cout<< " not found in PGA with filename: " + inFileName<< endl;
    }
    counter = 0;
    //Start reading in values.
    f_in >> time >> height;
    while(f_in) {
        if (counter < wSize) {
            wave.push_back(height); // Same wave
            counter++;
        }
        else {
            basel = 0.0;
            //Calculate and subtract the baseline (For LUNA results this looks like around 2244?).
            for (int i = 0; i < baseLEnd; i++) {
                basel += wave[i]/baseLEnd;
            }
            for (int i = 0; i < wave.size(); ++i) {
                wave[i] -= basel;
            }
            //Find amplitude value and sample value for the wave.
            amplitudeVal = maxModVal(wave);
            sampleVal = wave[sampleNo];
            PGAVal = abs(sampleVal - amplitudeVal);
            //print to output
            if (f_out.is_open()) {
                f_out << PGAVal << endl;
            } else {
                cout << "Unable to open file: " + outFileName << endl;
            }
            wave.clear();
            counter = 0;
        }
        f_in >> time >> height;
    }
    f_in.close();
    f_out.close();
    cout<<"                       PGA Completed                    "<<endl;

}
//-------------------------------------------------Run Comparison Methods-----------------------------------------------
//It's become necessary to compare various aspects of runs to determine what is causing the gradual increase in
//neutron rates with real time. The earliest runs in real time from LUNA are dump_001_wf_0 and dump_001_wf_1.
//----------------------------------------------------------------------------------------------------------------------

//Method to calculate the average peak height for the waves for 2 runs for comparison. The baseline is first subtracted
//from the wave, then the peak value is extracted
void peakValAverage(string inFileName, string outFileName, int wSize, int baseLEnd){
    int counter, time;
    vector<double> peakHeights, wave;
    double height, avgHeight, basel;
    fstream f_in;
    f_in.open(inFileName.c_str(),std::fstream::in);
    if(!f_in){
        cout<< " not found in peakValAverage with filename: " + inFileName << endl;
    }
    avgHeight = 0;
    counter = 0;
    //Start reading in values.
    f_in >> time >> height;
    while(f_in){
        if(counter<wSize){
            wave.push_back(height); // Same wave
            counter++;
        }
        else {
            basel = 0.0;
            //Find and subtract the baseline for the wave
            for(int i=0;i<baseLEnd;i++){
                basel+=wave[i]/baseLEnd;
            }
            for(int i=0;i<wSize;i++){
                wave[i]-=basel;
            }
            //Find maxVal for the wave.
            peakHeights.push_back(modMaxModVal(wave));
            wave.clear();
            counter = 0;
        }
        f_in >> time >> height;
    }
    f_in.close();
    for(int i=0;i<peakHeights.size();++i){
        avgHeight+=peakHeights[i]/peakHeights.size();
    }

    ofstream f_out(outFileName, ios::out | ios::app);
    if (f_out.is_open()) {
        f_out << inFileName<<" "<< avgHeight << endl;
    } else {
        cout << "Unable to open file: " + outFileName << endl;
    }
    f_out.close();
    cout<<"                       peakValAverage Completed                    "<<endl;

}


//Method to compare the average baseline for 2 runs for comparison.
void baselineAverage(string inFileName, string outFileName, int wSize, int baseLEnd){

    int counter, time;
    vector<double> baseLVals, wave;
    double height, avgBaseL, basel;
    fstream f_in;
    f_in.open(inFileName.c_str(),std::fstream::in);
    if(!f_in){
        cout<< " not found in baselineAverage with filename: " + inFileName << endl;
    }
    avgBaseL = 0;
    counter = 0;
    //Start reading in values.
    f_in >> time >> height;
    while(f_in){
        if(counter<wSize){
            wave.push_back(height); // Same wave
            counter++;
        }
        else {
            basel = 0.0;
            //Find the baseline for the wave
            for(int i=0;i<baseLEnd;i++){
                basel+=wave[i]/baseLEnd;
            }
            //Find maxVal for the wave.
            baseLVals.push_back(basel);
            wave.clear();
            counter = 0;
        }
        f_in >> time >> height;
    }
    f_in.close();
    for(int i=0;i<baseLVals.size();++i){
        avgBaseL+=baseLVals[i]/baseLVals.size();
    }

    ofstream f_out(outFileName, ios::out | ios::app);
    if (f_out.is_open()) {
        f_out << inFileName<<" "<< avgBaseL << endl;
    } else {
        cout << "Unable to open file: " + outFileName << endl;
    }
    f_out.close();
    cout<<"                       baselineAverage Completed                    "<<endl;
}

//Method to compare the average width for the neutron region and the low non-neutron region for 2 runs.
void regionWidthComparison(string inFileName1, string inFileName2, double lowThreshold, double highThreshold){
    vector<double> neutronVec1, nonNeutronVec1, neutronVec2, nonNeutronVec2;
    double width1, width2, avgNeutWidth1, avgNonWidth1, avgNeutWidth2, avgNonWidth2;
    avgNeutWidth1 = avgNeutWidth2 = avgNonWidth1 = avgNonWidth2 = 0.0;
    //First run
    fstream f_in1;
    f_in1.open(inFileName1.c_str(),std::fstream::in);
    if(!f_in1){
        cout<< " not found in regionWidthComparison first filename with filename: " + inFileName1 << endl;
    }
    //Start reading in values.
    f_in1 >> width1;
    while(f_in1){
        if((width1>lowThreshold)&&(width1<highThreshold)){
            neutronVec1.push_back(width1);
        }else if(width1<lowThreshold){
            nonNeutronVec1.push_back(width1);
        }
        f_in1 >> width1;
    }
    f_in1.close();

    for (int i=0;i<neutronVec1.size();++i){
        avgNeutWidth1+=neutronVec1[i];
    }
    avgNeutWidth1/=neutronVec1.size();

    for (int i=0;i<nonNeutronVec1.size();++i){
        avgNonWidth1+=nonNeutronVec1[i];
    }
    avgNonWidth1/=nonNeutronVec1.size();

    //Repeat for other run.
    fstream f_in2;
    f_in2.open(inFileName2.c_str(),std::fstream::in);
    if(!f_in2){
        cout<< " not found in regionWidthComparison second filename with filename: " + inFileName2 << endl;
    }
    //Start reading in values.
    f_in2 >> width2;
    while(f_in2){
        if((width2>lowThreshold)&&(width2<highThreshold)){
            neutronVec2.push_back(width2);
        }else if(width2<lowThreshold){
            nonNeutronVec2.push_back(width2);
        }
        f_in2 >> width2;
    }
    f_in2.close();
    //Average the values over their vectors.
    for (int i=0;i<neutronVec2.size();++i){
        avgNeutWidth2+=neutronVec2[i];
    }
    avgNeutWidth2/=neutronVec2.size();

    for (int i=0;i<nonNeutronVec2.size();++i){
        avgNonWidth2+=nonNeutronVec2[i];
    }
    avgNonWidth2/=nonNeutronVec2.size();

    //Print em.
    cout<<"For the input file "<<inFileName1<<" the average neutron region width is: "<<avgNeutWidth1<<endl;
    cout<<"and the average low region non-neutron width is: "<<avgNonWidth1<<endl;
    cout<<"For the input file "<<inFileName2<<" the average neutron region width is: "<<avgNeutWidth2<<endl;
    cout<<"and the average low region non-neutron width is: "<<avgNonWidth2<<endl<<endl;
    cout<<"                       regionWidthComparison Completed                    "<<endl;

}

//method to calculate the average deviation from the baseline for diagnosing electronic noise in LUNA runs.
void baselineDeviation(string inFileName, string outFileName, int wSize, int baseLEnd){
    fstream f_in;
    f_in.open(inFileName.c_str(),std::fstream::in);
    if(!f_in){
        cout<< " not found in baselineDeviation with filename: " + inFileName << endl;
    }
    vector<double> wave;
    double height, basel, deviation, delta;
    int counter = 0, time;
    f_in >> time >> height;
    while(f_in){
        if(counter<wSize){
            wave.push_back(height); // Same wave
            counter++;
        }
        else {
            basel = 0.0;
            //Find the baseline for the wave
            for(int i=0;i<baseLEnd;i++){
                basel+=wave[i]/baseLEnd;
            }
            for(int i=0;i<baseLEnd;i++){
                delta = wave[i]-basel;
                deviation+=delta*delta/baseLEnd;
            }
            deviation = sqrt(deviation);
            wave.clear();
            counter = 0;
        }
        f_in >> time >> height;
    }
    ofstream f_out(outFileName, ios::out | ios::app);
    //cout << inFileName <<" "<<deviation<< endl;
    if (f_out.is_open()) {
        f_out << inFileName <<" "<<deviation<< endl;
    } else {
        cout << "Unable to open file: " + outFileName << endl;
    }
    f_out.close();
    f_in.close();;
    cout<<"                       baselineDeviation Completed                    "<<endl;
}

//Method to sort the LUNA runs by detector.
void sortedLUNA(string inFileName, string outFileName0, string outFileName1){
    fstream f_in;
    f_in.open(inFileName.c_str(),std::fstream::in);
    if(!f_in){
        cout<< " not found in sortedLUNA with filename: " + inFileName << endl;
    }
    string filename1, filename0;
    double LUNAVal1, LUNAVal0;
    //Start reading in values.
    ofstream f_out0(outFileName0, ios::out | ios::app), f_out1(outFileName1, ios::out | ios::app);
    f_in >> filename0 >> LUNAVal0 >>filename1>>LUNAVal1;
    while(f_in){
        if (f_out0.is_open()) {
            f_out0 << filename0 <<" "<<LUNAVal0<<endl;
        } else {
            cout << "Unable to open file0 " << endl;
        }
        if (f_out1.is_open()) {
            f_out1 << filename1 <<" "<<LUNAVal1<<endl;
        } else {
            cout << "Unable to open file1 " << endl;
        }
        f_in >> filename0 >> LUNAVal0 >>filename1>>LUNAVal1;
    }
    f_out0.close();
    f_out1.close();
    f_in.close();
    cout<<"                       sortedLUNA Completed                    "<<endl;
}


//------------------------------------------Derived Quantities----------------------------------------------------------
//Section to calculate the physically derived quantities of the detectors/neutrons coming through from the
//above methods. Quantites such as flux, efficiency etc.
//----------------------------------------------------------------------------------------------------------------------
//Method to sort the neutron cases from the non-neutrons. Outputs the number of neutrons, non-neutrons,
//their rates per hour and per second, total events, the neutron flux and its error in cm^-2s^-1 and the
//difference between the rates of non-neutrons and neutrons with its error in s^-1. All derived from the
//Widths methods above.

void printWidthsDerivedQuantities(string inFileName, double lowThreshold, double highThreshold, double time){
    cout<<inFileName<<endl;
    fstream f_in;
    f_in.open(inFileName.c_str(),std::fstream::in);
    if(!f_in){
        cout<< " not found in printWidthsDerivedQuantities with filename: " + inFileName << endl;
    }
    double inVal, flux, fluxError;
    int numRejections = 0;
    int numNeutrons = 0;
    int total = 0;
    //Start reading in values.
    f_in >> inVal;
    while(f_in){
        if((!(inVal<lowThreshold))&&(!(inVal>highThreshold))){
            numNeutrons++;
        }
        else{
            numRejections++;
        }
        total++;
        f_in>>inVal;
    }

    //calculate errors and other associated values
    double A = EJ426DETY*EJ426DETX;
    double AT = EJ426DETY*EJ426DETX*time;
    flux = numNeutrons/AT;
    fluxError = sqrt(numNeutrons/(AT*AT)
                     + numNeutrons*numNeutrons*DETAREAERROR*DETAREAERROR/(AT*AT*A*A)
                     + numNeutrons*numNeutrons*TIMEERR*TIMEERR/(AT*AT*time*time));
    double neutRateErr = sqrt(numNeutrons*(1+numNeutrons*TIMEERR*TIMEERR/(time*time)))/time;
    double nonNeutRateErr = sqrt(numRejections*(1+numRejections*TIMEERR*TIMEERR/(time*time)))/time;
    double neutRateHrErr = sqrt(numNeutrons*(1+numNeutrons*TIMEERR*TIMEERR/(time*time)))*3600/time;
    double nonNeutRateHrErr = sqrt(numRejections*(1+numRejections*TIMEERR*TIMEERR/(time*time)))*3600/time;

    cout<<"For the input run widths file: "<<inFileName<<endl
    <<"which was "<<time<<"s long"<<endl
    <<"using the low threshhold value of "<<lowThreshold<< " and a high of "<< highThreshold
    <<" to discriminate widths, the number of neutrons is: "<<numNeutrons<<endl
    <<"and the number of rejections is: "<<numRejections<<endl
    <<"with a total number of events: "<<total<<endl
    <<"this results in a calculated neutron flux, for the detector, of: "<<flux<<"cm^-2s^-1"<<endl
    <<"with an associated error of: "<<fluxError<<"cm^-2s^-1"<<endl
    <<"the neutron rate for this run is: "<<numNeutrons/time<<"s^-1 with an error of: "<<neutRateErr<<"s^-1" <<endl
    <<"which, in units of hours, is: "<<numNeutrons*3600/time<<"hr^-1 with an error of: "<<neutRateHrErr<<"hr^-1" <<endl
    <<"the non-neutron rate for this run is: "<<numRejections/time<<"s^-1 with an error of: "<<nonNeutRateErr<<"s^-1" <<endl
    <<"which, in units of hours, is: "<<numRejections*3600/time<<"hr^-1 with an error of: "<<nonNeutRateHrErr<<"hr^-1" <<endl
    <<"The difference between non-neutron and neutron rates is: "<<numRejections/time - numNeutrons/time<<"s^-1 with "
    <<"error: "<<sqrt(neutRateErr*neutRateErr + nonNeutRateErr*nonNeutRateErr)<<"s^-1" <<endl;
    f_in.close();
    cout<<"                       printWidthsDerivedQuantities Completed                    "<<endl;

}

void printWidthsDerivedQuantitiesOutFile(string inFileName, string outFileName, double lowThreshold, double highThreshold, double time){
    cout<<inFileName<<endl;
    fstream f_in;
    f_in.open(inFileName.c_str(),std::fstream::in);
    if(!f_in){
        cout<< " not found in printWidthsDerivedQuantitiesOutFile with filename: " + inFileName << endl;
    }
    double inVal;
    int numRejections = 0;
    int numNeutrons = 0;
    int total = 0;
    //Start reading in values.
    f_in >> inVal;
    while(f_in){
        if((!(inVal<lowThreshold))&&(!(inVal>highThreshold))){
            numNeutrons++;
        }
        else{
            numRejections++;
        }
        total++;
        f_in>>inVal;
    }
    ofstream f_out(outFileName, ios::out | ios::app);
    if (f_out.is_open()) {
        f_out << inFileName << " "<< time<<" "<<numNeutrons<<endl;
    } else {
        cout << "Unable to open file: " + outFileName << endl;
    }
    f_out.close();
    f_in.close();
    cout<<"                       printWidthsDerivedQuantitiesOutFile Completed                    "<<endl;

}


//Method to print out to a file the filename and the neutron rate for that run based on the thresholds given. The thresholds
//are to be decided upon inspection of the widths files produced by the above methods.
void WidthDerivedNeutronRate(string inFileName, string outFileName, double lowThreshold, double highThreshold, double time){
    cout<<inFileName<<endl;
    fstream f_in;
    f_in.open(inFileName.c_str(),std::fstream::in);
    if(!f_in){
        cout<< " not found in WidthDerivedNeutronRate with filename: " + inFileName << endl;
    }
    double inVal;
    int numNeutrons = 0;
    //Start reading in values.
    ofstream f_out(outFileName, ios::out | ios::app);
    f_in >> inVal;
    while(f_in){
        if((!(inVal<lowThreshold))&&(!(inVal>highThreshold))){
            numNeutrons++;
        }
        f_in>>inVal;
    }
    if (f_out.is_open()) {
        f_out << inFileName <<" "<<numNeutrons/time<<endl;
    } else {
        cout << "Unable to open file " << endl;
    }
    f_out.close();
    f_in.close();
    cout<<"For the input run widths file "<<inFileName<<" the neutron rate is: "<<numNeutrons/time<<endl;
    cout<<"                       WidthDerivedNeutronRate Completed                    "<<endl;
}

double widthDerivedNeutronRateVal(string inFileName, double lowThreshold, double highThreshold, double time){
    cout<<inFileName<<endl;
    fstream f_in;
    f_in.open(inFileName.c_str(),std::fstream::in);
    if(!f_in){
        cout<< " not found in WidthDerivedNeutronRateVal with filename: " + inFileName << endl;
    }
    double inVal;
    int numNeutrons = 0;
    //Start reading in values.
    f_in >> inVal;
    while(f_in){
        if((!(inVal<lowThreshold))&&(!(inVal>highThreshold))){
            numNeutrons++;
        }
        f_in>>inVal;
    }
    f_in.close();
    cout<<"                       WidthDerivedNeutronRateVal Completed                    "<<endl;
    return numNeutrons/time;
}

//Method to calculate the efficiency of the detector from the number of neutrons measured by the widths method (takes in a _Widths file),
// an orientation and the activity of the source. Will produce both the absolute and intrinsic effeciency. The input orientation must be
// either "horizontal" or "vertical", being the largest faces of the detector facing up and down or left and right respectively.
void WidthsDerivedEfficiencies(string inFileName, string outFileName, string orientation, double distanceInMetres, double sourceActivity,
                                double lowThreshold, double highThreshold, double time){
    if (!(orientation == "horizontal")&&!(orientation == "vertical")){
        cout<<"Please enter an orientation of 'vertical' or 'horizontal'."<<endl;
        exit(1);
    }
    fstream f_in;
    f_in.open(inFileName.c_str(),std::fstream::in);
    if(!f_in){
        cout<< " not found in WidthsDerivedEfficiency with filename: " + inFileName << endl;
    }
    double neutronRate = widthDerivedNeutronRateVal(inFileName, lowThreshold, highThreshold, time);
    double absoluteEfficiency = neutronRate/sourceActivity;
    double solidAngle, detectorWidth, detectorDepth, bonusDistance, trueDistance;
    if(orientation == "horizontal"){
        detectorWidth = DETZ/100; //metres
        detectorDepth = DETY/100; //metres
    }else{
        detectorWidth = DETY/100;
        detectorDepth = DETZ/100;
    }
    bonusDistance = detectorDepth/2;
    trueDistance = distanceInMetres + bonusDistance;
    solidAngle = 4*atan(detectorWidth*EJ426DETX/
                                (4*trueDistance*sqrt(detectorWidth*detectorWidth/
                                                             4+EJ426DETX*EJ426DETX/4+trueDistance*trueDistance)));
    double intrinsicEfficiency = absoluteEfficiency*4*M_PI/solidAngle;

    ofstream f_out(outFileName, ios::out | ios::app);
    if (f_out.is_open()) {
        f_out << inFileName <<" "<<absoluteEfficiency<<" "<< intrinsicEfficiency<< endl;
    } else {
        cout << "Unable to open file: " + outFileName << endl;
    }
    f_out.close();
    f_in.close();;
    cout<<"                       WidthsDerivedEfficiency Completed                    "<<endl;
}


//----------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------MAIN-----------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------
int main() {

    reprint("AmBe_Spectrum.txt","AmBe_Spectrum_Processed.txt");
    int wSize, //Number of points in the waveform.
            baseLEnd, //Point up to which only the baseline is present.
            tailW, //Point at which the tail of the pulse ends.
            peakXValue, //Point at which the peak value of the wave is.
            wStart, //Point at which the wave starts, in practice is the same as baseLEnd.
            wEnd, //Point at which the pulse ends, in practice is the same as tailW.
            PGASampleVal, //Sample value for the PGA method.
            widthLowCut, //Low cut point for the method counting the number of neutrons.
            widthHighCut; //High cut point for the method counting the number of neutrons.
    double sourceDistance, //Distance from the detector to the source, if there is one, for background runs this should
                           //be set to 0 in the input file.
            AmBeSourceActivity, //Neutron source activity.
            runTime; //Duration of the run in seconds.
    string fileModifier, //Type of input file used (.txt, .dat, .csv etc.)
            orientation, //Horizontal or vertical detector orientation?
            location, //Location of the detector runs that sets up other variables
            filename, //Run name.
            fileDestination; //Folder containing the files.
    string fileDetails = "File Details.txt"; //Name of the input file containing all details of the runs.
    fstream f_in;
    f_in.open(fileDetails.c_str(),std::fstream::in);
    if(!f_in){
        cout<< " not found MAIN" << endl;
    }
    //Start reading in values. They are all stored in an input file to be edited upon reception of new data.
    f_in >> filename >> runTime >> location >> fileDestination >> sourceDistance >> orientation;
    int counter = 0;
    cout<<endl<<endl<<"        ---------------Beginning PSD Codes--------------- "<<endl;
    while(f_in){
        counter ++;
        if ((counter<4)||(counter>111)){
            f_in  >> filename >> runTime >> location >> fileDestination >> sourceDistance >> orientation;
            continue;
        }

        cout<<"             Starting at line "<<counter<<" in "<<fileDetails<<endl;
        if(location == "SeptEdinburgh"){
            wSize = 1000;
            baseLEnd = 100;
            tailW = 600;
            peakXValue = 200;
            wStart = 100;
            wEnd = 800;
            PGASampleVal = 600;
            fileModifier = ".csv";
            AmBeSourceActivity = 2.738E5; //Neutrons per second
            widthLowCut = 5;
            widthHighCut = 50;
        }else if(location == "LUNA"){
            wSize = 4000;
            baseLEnd = 30;
            tailW = 200;
            peakXValue = 34;
            wStart = 30;
            wEnd = 100;
            PGASampleVal = 100;
            fileModifier = ".dat";
        }else if(location == "JanEdinburgh"){
            wSize = 100000;
            baseLEnd = 10000;
            tailW = 38000;
            peakXValue = 22000;
            wStart = 19000;
            wEnd = 60000;
            PGASampleVal = 40000;
            fileModifier = ".txt";
            AmBeSourceActivity = 2.737E5; //Neutrons per second
            widthLowCut = 19000;
            widthHighCut = 40000;
        }else if(location == "FebEdinburgh"){
            wSize = 10000;
            baseLEnd = 1000;
            tailW = 3800;
            peakXValue = 2200;
            wStart = 1900;
            wEnd = 6000;
            PGASampleVal = 4000;
            fileModifier = ".txt";
            AmBeSourceActivity = 2.737E5; //Neutrons per second
            widthLowCut = 1900;
            widthHighCut = 4000;
        }else{
            cout<< "Please make sure the file format contains \"LUNA\", \"JanEdinburgh\" or \"FebEdinburgh\""<<endl;
            break;
        }
        cout << "Filename: "<< filename << ", runTime: "<< runTime << "s, location: "
        << location << ", fileDestination: " << fileDestination << " "<<endl
        << "sourceDistance: "<< sourceDistance << "m, orientation: " << orientation << endl;

        Widths(fileDestination + filename + fileModifier,
               fileDestination + "Widths/" + filename + "_Widths.txt", 0.5, wSize, baseLEnd);

        printWidthsDerivedQuantities(fileDestination + "Widths/" + filename + "_Widths.txt", widthLowCut,
                                     widthHighCut, runTime);

        printWidthsDerivedQuantitiesOutFile(fileDestination + "Widths/" + filename + "_Widths.txt",
                                            fileDestination + "Derived Quantities/timesandnumneutrons.txt",
                                            widthLowCut, widthHighCut, runTime);

        numWaves(fileDestination + filename + fileModifier, wSize);

        widthBinTimeNormalised(fileDestination + "Widths/" +filename + "_Widths.txt",
                               fileDestination + "Time Normalised/time_normalised_" + filename + "_Widths.txt",
                               runTime, 1, wSize);

        totalIntVsWidth(fileDestination + filename + fileModifier, fileDestination + "Total Integral vs Width/" +
                filename + "_Total_Integral_vs_Widths.txt", 0.5, wSize, baseLEnd, wStart, wEnd );

        peakTailIntegrate(fileDestination + filename + fileModifier, fileDestination + "Tail vs Peak Integral/"
                                                                     + filename + "_Tail_vs_Peak_Integral.txt",wSize,
                                                                     baseLEnd, peakXValue, tailW);

        PGA(fileDestination + filename + fileModifier, fileDestination + "PGA/" + filename + "_PGA.txt",
            PGASampleVal, wSize, baseLEnd);

        firstTen(fileDestination + filename + fileModifier,
                 fileDestination + "First Ten/" + filename + "_First Ten.txt", wSize);

        baselineAdjust(fileDestination + filename + fileModifier,
                       fileDestination + "Baseline Adjusted/" + filename + "_Baseline Adjusted.txt", wSize, baseLEnd);

        totalIntVsWidthPostBaselineAdjusted(fileDestination + "Baseline Adjusted/" + filename + "_Baseline Adjusted.txt",
                                                fileDestination + "Total Integral vs Width PBLA/" + filename +
                                                "_Total_Integral_vs_Width.txt", 0.5, wSize, wStart, wEnd);

        if((location=="SeptEdinburgh")||(location=="JanEdinburgh")||(location=="FebEdinburgh")){
            WidthDerivedNeutronRate(fileDestination + "Widths/" + filename + "_Widths.txt",
                                    fileDestination + "Derived Quantities/FWHM_derived_neutron_rate.txt",
                                    widthLowCut, widthHighCut, runTime);

            WidthsDerivedEfficiencies(fileDestination + "Widths/" + filename + "_Widths.txt",
                                      fileDestination + "Derived Quantities/FWHM_derived_neutron_absolute_and_intrinsic_efficiency.txt",
                                      orientation, sourceDistance, AmBeSourceActivity, widthLowCut, widthHighCut, runTime);

        }

        baselineDeviation(fileDestination + filename + fileModifier,
                          fileDestination + "Derived Quantities/Baseline Deviation.txt", wSize, baseLEnd);

        //peakValAverage(fileDestination+filename+fileModifier,fileDestination+"Derived Quantities/AvgPeak.txt", wSize, baseLEnd);
        baselineAverage(fileDestination+filename+fileModifier,fileDestination+"Derived Quantities/AvgBasel.txt", wSize, baseLEnd);
        cout << endl;
        f_in  >> filename >> runTime >> location >> fileDestination >> sourceDistance >> orientation;
    }
    f_in.close();

    //sortedLUNA("LUNA/Derived Quantities/Baseline Deviation.txt", "LUNA/Derived Quantities/Baseline Deviation 0.txt",
    //                 "LUNA/Derived Quantities/Baseline Deviation 1.txt");
    sortedLUNA("LUNA/Derived Quantities/AvgPeak.txt", "LUNA/Derived Quantities/AvgPeak0.txt",
               "LUNA/Derived Quantities/AvgPeak1.txt");
    sortedLUNA("LUNA/Derived Quantities/AvgBasel.txt", "LUNA/Derived Quantities/AvgBasel0.txt",
               "LUNA/Derived Quantities/AvgBasel1.txt");

    FoM(75, 2, 5, 2, 50, 2);
    //Need to make some comparison between runs 1 and 7 in air for both detectors. They are exhibiting
    //a steady increase in neutron count and a flat progression in non-neutron count. This could potentially
    //be radon, but other investigations must be performed to rule out certain options.

    peakValComparison("LUNA/dump_001_wf_0.dat", "LUNA/dump_007_wf_0.dat", wSize, baseLEnd);
    peakValComparison("LUNA/dump_001_wf_1.dat", "LUNA/dump_007_wf_1.dat", wSize, baseLEnd);

    baselineComparison("LUNA/dump_001_wf_0.dat", "LUNA/dump_007_wf_0.dat", wSize, baseLEnd);
    baselineComparison("LUNA/dump_001_wf_1.dat", "LUNA/dump_007_wf_1.dat", wSize, baseLEnd);

    regionWidthComparison("LUNA/dump_001_wf_0_Widths.dat", "LUNA/dump_007_wf_0_Widths.dat", 7.0, 50.0);
    regionWidthComparison("LUNA/dump_001_wf_1_Widths.dat", "LUNA/dump_007_wf_1_Widths.dat", 7.0, 50.0);

    regionWidthComparison("LUNA/poly_000_wf_0_Widths.dat", "LUNA/poly_008_wf_0_Widths.dat", 7.0, 50.0);
    regionWidthComparison("LUNA/poly_000_wf_1_Widths.dat", "LUNA/poly_008_wf_1_Widths.dat", 7.0, 50.0);

    regionWidthComparison("LUNA/AmBe_002_wf_0_Widths.dat", "LUNA/AmBe_010_wf_0_Widths.dat", 7.0, 50.0);
    regionWidthComparison("LUNA/AmBe_002_wf_1_Widths.dat", "LUNA/AmBe_010_wf_1_Widths.dat", 7.0, 50.0);

    printf("\a");
    return 0;
}