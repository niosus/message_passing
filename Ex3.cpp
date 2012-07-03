#include <NMath.h>
#include <vector>
#include <CVector.h>
#include <string>

using namespace std;
using namespace NMath;

#define bigNumber 99999999


void vectOut(vector<float> out)
{
	for (int i=0; i<out.size();++i)
	{
		cout<<out[i]<<" ";
	}
	cout<<endl;
}

void vectOut(vector<int> out)
{
	for (int i=0; i<out.size();++i)
	{
		cout<<out[i]<<" ";
	}
	cout<<endl;
}

vector<vector<float> > CreateUnaryCosts(
		vector<float>& pixels1,
		vector<float>& pixels2,
		int L)
{
	vector<vector<float> > result(pixels1.size());
	if (pixels1.size()!=pixels2.size())
	{
		string error="error: wrong dimensions in CreateUnaryCosts";
		throw error;
		return result;
	}
	for (int x=0; x<pixels1.size(); ++x)
	{
		vector<float> costs(L);
		for (int disp=0; disp<L; ++disp)
		{
			if (x-disp>=0)
				costs[disp]=abs(pixels1[x]-pixels2[x-disp]);
			else
				costs[disp]=bigNumber;
		}
		result[x]=costs;
	}
	return result;
}


vector<float> CreateMessage(
		const vector<float>& unaryCosts,
		const vector<float>& prevMessage,
		float binaryWeight)
{
	bool prevMessageExists=false;
	if (prevMessage.size()==unaryCosts.size())
	{
		prevMessageExists=true;
	}
	int L=unaryCosts.size(); //number of labels
	float tempCost=0;
	float binaryCost=0;
	vector<float> message(L);
	for (int j=0; j<L; ++j)
	{
		float minCost=bigNumber;
		for (int i=0; i<L; ++i)
		{
			if (i==j)
			{
				binaryCost=0;
			}
			else
			{
				binaryCost=binaryWeight; //binary weight corresponds to Potts model, int (float)
			}
			if (prevMessageExists)
				tempCost=unaryCosts[i]+prevMessage[i]+binaryCost;
			else
				tempCost=unaryCosts[i]+binaryCost;
			minCost=min(minCost,tempCost);
		}
		message[j]=minCost;
	}
	return message;
}

vector<vector<float> > CreateMessageArray(vector<float>& pixels1, vector<float>& pixels2, int L, int direction)
{
	vector<vector<float> > messages(pixels1.size());
	vector<vector<float> > unaryCosts=CreateUnaryCosts(pixels1,pixels2,L);
	vector<float> message;
	float binaryWeight=40;
	if (direction>0)
	{
		for (int i=0; i<pixels1.size(); ++i)
		{
			message=CreateMessage(unaryCosts[i],message,binaryWeight);
			messages[i]=message;
		}
	}
	else
	{
		for (int i=pixels1.size()-1; i>=0; --i)
		{
			message=CreateMessage(unaryCosts[i],message,binaryWeight);
			messages[i]=message;
		}
	}
	cout<<messages.size()<<endl;
	return messages; //an array of messages for each node
}

vector<int> DecideLabels(
		const vector<vector<float> >& mesDir1,
		const vector<vector<float> >& mesDir2,
		const vector<vector<float> >& unaryCosts)
{
	int numOfPixels=mesDir1.size();
	//cout<<numOfPixels<<endl;
	vector<int> result(numOfPixels);
	for (int i=0; i<numOfPixels; ++i)
	{
		float minVal=bigNumber;
		for (int j=0; j<mesDir1[i].size(); ++j)
		{
			float temp=mesDir1[i][j]+mesDir2[i][j]+unaryCosts[i][j];
			if (minVal>temp)
			{
				minVal=temp;
				result[i]=j;
			}
		}
	}
	return result;
}

vector<float> CVect2Vect(const CVector<float>& cvect)
{
	vector<float> result;
	//cout<<cvect.size()<<endl;
	for (int i=0; i<cvect.size(); ++i)
	{
		result.push_back(cvect(i));
	}
	return result;
}



int main()
{
	CMatrix<float> imageL;
	imageL.readFromPGM("tsukubaL.pgm");
	CMatrix<float> imageR;
	imageR.readFromPGM("tsukubaR.pgm");
	CMatrix<float> result(imageL.xSize(),imageL.ySize());
	vector<vector<float> > MesDir1, MesDir2, UnaryCosts;
	CVector<float> CpixelsL, CpixelsR;
	vector<float> pixelsL, pixelsR;
	vector<int> temp_res;
	for (int y=0; y<imageL.ySize(); ++y)
	{
		cout<<"row = "<< y << endl;
		for (int x=0; x<imageL.xSize(); ++x)
		{
			pixelsL.push_back(imageL(x,y));
			pixelsR.push_back(imageR(x,y));
		}

		int maxDisparity=15;
		int direction=1;
		MesDir1=CreateMessageArray(pixelsL,pixelsR,maxDisparity,direction);
		direction=-1;
		MesDir2=CreateMessageArray(pixelsL,pixelsR,maxDisparity,direction);
		UnaryCosts=CreateUnaryCosts(pixelsL,pixelsR,maxDisparity);
		temp_res=DecideLabels(MesDir1,MesDir2,UnaryCosts);
		vectOut(temp_res);
		for (int x=0; x<imageL.xSize(); ++x)
		{
			result(x,y)=temp_res[x];
		}
		pixelsL.clear();
		pixelsR.clear();
	}
	result.normalize(0,255,0,15);
	result.writeToPGM("result.pgm");
	return 0;
}
