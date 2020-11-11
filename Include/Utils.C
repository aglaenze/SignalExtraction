#include <iostream>
#include <fstream>
#include <dirent.h>
#include <string>
#include <cmath>


using namespace std;

bool Initiate(string period, Double_t& mMin, Double_t& mMax, Double_t& ptMin, Double_t& ptMax, bool& useCuts, bool& logScale, bool& drawPulls, bool& exp, bool& exclusiveOnly) {
	ifstream file("input-" + period + ".txt", ios::in);
	string a, b;
	bool mMinFound = false, mMaxFound = false, ptMinFound = false, ptMaxFound = false, useCutsFound = false, logScaleFound = false, drawPullsFound = false, expFound = false, exclusiveOnlyFound = false;
	if (file) {
		string line {};
		getline(file, line);	//first line does not contains info
		cout << "#################################" << endl;
		while(getline(file, line)) {
			if (line.find("mMin") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> mMin;
				mMinFound = true;
				cout << "#\tmMin = " << mMin << "\t\t#" << endl;
			}
			else if (line.find("mMax") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> mMax;
				mMaxFound = true;
				cout << "#\tmMax = " << mMax << "\t\t#" << endl;
			}
			else if (line.find("ptMin") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> ptMin;
				ptMinFound = true;
				cout << "#\tptMin = " << ptMin << "\t\t#" << endl;
			}
			else if (line.find("ptMax") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> ptMax;
				ptMaxFound = true;
				cout << "#\tptMax = " << ptMax << "\t\t#" << endl;
			}
			else if (line.find("useCuts") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> useCuts;
				useCutsFound = true;
				cout << "#\tuseCuts = " << useCuts << "\t\t#" << endl;
			}
			else if (line.find("logScale") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> logScale;
				logScaleFound = true;
				cout << "#\tlogScale = " << logScale << "\t\t#" << endl;
			}
			else if (line.find("drawPulls") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> drawPulls;
				drawPullsFound = true;
				cout << "#\tdrawPulls = " << drawPulls << "\t\t#" << endl;
			}
			else if (line.find("exp") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> exp;
				expFound = true;
				cout << "#\texp = " << exp << "\t\t\t#" << endl;
			}
			else if (line.find("exclusiveOnly") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> exclusiveOnly;
				exclusiveOnlyFound = true;
				cout << "#\texclusiveOnly = " << exclusiveOnly << "\t#" << endl;
			}
		}
		cout << "#################################" << endl;
		cout << endl;
		return (mMinFound && mMaxFound && ptMinFound && ptMaxFound && useCutsFound && logScaleFound && drawPullsFound && expFound && exclusiveOnlyFound);
	}
	else cout << "Error: not possible to open input.txt file in reading mode" << endl;
	return false;
}


bool GetParameters(string period, Double_t& bExc, Double_t& gammaPbYield) {
	ifstream file("input-" + period + ".txt", ios::in);
	string a, b;
	bool bExcFound = false, gammaPbYieldFound = false;
	if (file) {
		string line {};
		getline(file, line);	//first line does not contains info
		while(getline(file, line)) {
			if (line.find("bExc") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> bExc;
				//cout << "\n\n\n\n\n\nb = " << bExc << endl << endl << endl << endl << endl;
				bExcFound = true;
			}
			else if (line.find("gammaPbYield") != string::npos) {
				stringstream stream(line);
				stream >> a >> b >> gammaPbYield;
				gammaPbYieldFound = true;
			}
		}
		return (bExcFound && gammaPbYieldFound);
	}
	else cout << "Error: not possible to open input.txt file in reading mode" << endl;
	return false;
}


