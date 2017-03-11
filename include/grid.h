#ifndef GRID_H
#define GRID_H
#include <fstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>
using namespace std;
class grid{
//Reads external grid files, computes additional grid data, holds geometric data
//Any particular flow problem should only need one grid
    typedef std::vector< std::vector<double> > vec2D;
    public:
        grid(std::string xfilename, std::string yfilename){
            loadGrid(xfilename, yfilename); //sets x and y corner points
            defineCenterPoints();
            defineSideLengths();
            defineSideNormals();
			defineIJNormals();
            defineAreas();
        }
        vec2D xCorner, yCorner, xCenter, yCenter;
        vec2D xSside, ySside, xWside, yWside;
        vec2D xSnorm, ySnorm, xWnorm, yWnorm;
        vec2D xInorm, yInorm, xJnorm, yJnorm; //the centered, averaged normals
        vec2D area;
        int N, M; //number of points in the xi, and eta directions
    private:
        double **dx, **dy; //length parameters of cells between corners. Indices match edge-point indices.
        void loadGrid(std::string xfilename, std::string yfilename){ //Define N and M in here; x_corners = new double[][]; more initializing stuff
            std::ifstream xfile, yfile;
            std::string readin;
            int lineNum=0;
            xfile.open(xfilename.c_str());
            yfile.open(yfilename.c_str());
            if(!xfile.good()){std::cout<<"Error Reading x-grid file!\n";}
            if(!yfile.good()){std::cout<<"Error Reading y-grid file!\n";}

            //Figure out size of the csv table
            getline(xfile, readin);
            lineNum++;
            M=std::count(readin.begin(),readin.end(),',')+1;
            while(xfile.good()){
                getline(xfile, readin);
                lineNum++;
            }
            N=lineNum-1;
            xfile.clear();
            xfile.seekg(0, xfile.beg);//restart xfile from the beginning

            //Start Loading in data now that we know how much information there is
            //Allocate the correct amount of memory for xCorners and yCorners
            xCorner.resize(N, std::vector<double>(M));
            yCorner.resize(N, std::vector<double>(M));

            for(int i=0; i<N; i++){
                //read all data entries except the last in the line (the last line is the only one that doesnt have a comma at the end of it)
                for(int j=0; j<M-1; j++){
                    getline(xfile,readin,',');
                    xCorner[i][j]=atof(readin.c_str());
                    getline(yfile,readin,',');
                    yCorner[i][j]=atof(readin.c_str());
                }
                getline(xfile,readin);
                xCorner[i][M-1]=atof(readin.c_str());
                getline(yfile,readin);
                yCorner[i][M-1]=atof(readin.c_str());
            }
            xfile.close();
            yfile.close();
        }
        void defineCenterPoints(){ //interpolate cell centers
            xCenter.resize(N-1, std::vector<double>(M-1));
            yCenter.resize(N-1, std::vector<double>(M-1));
            for(int i=0; i<N-1; i++){
                for(int j=0; j<M-1; j++){
                    xCenter[i][j]=(xCorner[i][j]+xCorner[i][j+1]+xCorner[i+1][j]+xCorner[i+1][j+1])/4;
                    yCenter[i][j]=(yCorner[i][j]+yCorner[i][j+1]+yCorner[i+1][j]+yCorner[i+1][j+1])/4;
                }
            }
        }
        void defineAreas(){
            area.resize(N-1, std::vector<double>(M-1));
            for(int i=0; i<N-1; i++){
                for(int j=0; j<M-1; j++){
                    area[i][j]=0.5*abs(xSside[i][j]*yWside[i][j]-xWside[i][j]*ySside[i][j]);
                }
            }
        }
        void defineSideLengths(){
            xSside.resize(N-1, std::vector<double>(M-1));
            ySside.resize(N-1, std::vector<double>(M-1));
            xWside.resize(N-1, std::vector<double>(M-1));
            yWside.resize(N-1, std::vector<double>(M-1));
            for(int i=0; i<N-1; i++){
                for(int j=0; j<M-1; j++){
                    xSside[i][j]=abs(xCorner[i+1][j]-xCorner[i][j]);
                    ySside[i][j]=abs(yCorner[i+1][j]-yCorner[i][j]);
                    xWside[i][j]=abs(xCorner[i][j+1]-xCorner[i][j]);
                    yWside[i][j]=abs(yCorner[i][j+1]-yCorner[i][j]);
                }
            }
        }
        void defineSideNormals(){
            xSnorm.resize(N-1, std::vector<double>(M-1));
            ySnorm.resize(N-1, std::vector<double>(M-1));
            xWnorm.resize(N-1, std::vector<double>(M-1));
            yWnorm.resize(N-1, std::vector<double>(M-1));
            for(int i=0; i<N-1; i++){
                for(int j=0; j<M-1; j++){
                    ySnorm[i][j]=-(xCorner[i+1][j]-xCorner[i][j]);
                    xSnorm[i][j]=(yCorner[i+1][j]-yCorner[i][j]);
                    yWnorm[i][j]=(xCorner[i][j+1]-xCorner[i][j]);
                    xWnorm[i][j]=-(yCorner[i][j+1]-yCorner[i][j]);
                }
            }
        }

        void defineIJNormals(){
            xInorm.resize(N-1, std::vector<double>(M-1));
            yInorm.resize(N-1, std::vector<double>(M-1));
            xJnorm.resize(N-1, std::vector<double>(M-1));
            yJnorm.resize(N-1, std::vector<double>(M-1));
            for(int i=0; i<N-2; i++){
                for(int j=0; j<M-2; j++){
                    xInorm[i][j]=0.5*(xWnorm[i][j]+xWnorm[i+1][j]);
                    yInorm[i][j]=0.5*(yWnorm[i][j]+yWnorm[i+1][j]);
                    xJnorm[i][j]=0.5*(xSnorm[i][j]+xSnorm[i][j+1]);
                    yJnorm[i][j]=0.5*(ySnorm[i][j]+ySnorm[i][j+1]);
                }
            }

			//To get the boundary values for these s's, Snorm must be found at the extreme adge
/*			for (int i = 0; i < N - 1; i++) {

				xInorm[i][M-1] = 0.5*(xWnorm[i][j] + xWnorm[i + 1][j]);
				yInorm[i][j] = 0.5*(yWnorm[i][j] + yWnorm[i + 1][j]);
				xJnorm[i][j] = 0.5*(xSnorm[i][j] + xSnorm[i][j + 1]);
				yJnorm[i][j] = 0.5*(ySnorm[i][j] + ySnorm[i][j + 1]);

			}

			for (int j = 0; j < M - 1; j++) {
		
				

				}
			}*/
        }
};

#endif // GRID_H
