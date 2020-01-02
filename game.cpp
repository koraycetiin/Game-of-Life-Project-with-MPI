/*
*Koray Cetin
*2016400081
*Compiling
*Working
*Periodic - Checkered
*/
#include <bits/stdc++.h> 
#include <mpi.h>
using namespace std;
int worldRank;
int worldSize;
int edgeSize;
const int NOT_RECEIVED = -1;
const int MASTER_PROCESS_RANK = 0;
vector<vector<vector<int>>> gameMaps;
vector<int> ups;
vector<int> downs;
vector<int> lefts;
vector<int> rights;
vector<int> corners;
int stage = 0;

int getWorldRank (string direction){
    //Get world rank of the process in given direction
    int processPerRow = (int)sqrt(worldSize);
    bool isOnTopRow = worldRank <= processPerRow;
    bool isOnLeftmostColumn = worldRank % processPerRow == 1;
    bool isOnRightmostColumn = worldRank % processPerRow == 0;
    bool isOnBottomRow = worldRank >= worldSize - processPerRow;
    bool isOnSECorner = worldRank == worldSize - 1;
    bool isOnSWCorner = worldRank == worldSize - processPerRow;
    bool isOnNECorner = worldRank == processPerRow;
    bool isOnNWCorner = worldRank == 1;
    if(direction == "N"){
        if(isOnTopRow) {
            return worldSize - processPerRow + worldRank - 1;
            //Process on the first row
        }
        return worldRank - processPerRow;
    } else if(direction == "NW"){
        if(isOnTopRow){
            if(isOnNWCorner) {
                return worldSize - 1;
                //Process on the SE corner
            }
            return worldSize - processPerRow + worldRank - 2;
        } else {
            if(isOnLeftmostColumn){
                return worldRank - 1;
                //Process on the rightmost column
            }
            return worldRank - processPerRow - 1;
        }
    } else if(direction == "NE"){
        if(isOnTopRow){
            if(isOnNECorner){
                return worldSize - processPerRow;
                //Process on the SW corner
            }
            return worldSize - processPerRow + worldRank;
        } else {
            if(isOnRightmostColumn) {
                return worldRank - 2 * processPerRow + 1;
                //Process on the leftmost column
            }
            return worldRank - processPerRow + 1;
        }
    } else if(direction == "W"){
        if(isOnLeftmostColumn) {
            return worldRank + processPerRow - 1;
            //Process on the rightmost column 
        }
        return worldRank - 1;
    } else if(direction == "E"){
        if(isOnRightmostColumn){
            return worldRank - processPerRow + 1;
            //Process on the leftmost column
        }
        return worldRank + 1;
    } else if(direction == "S"){
        if(isOnBottomRow){
            if(isOnSECorner){
                return processPerRow;
                //Process on the NE corner
            }
            return worldRank % processPerRow;
        } 
        return worldRank + processPerRow;
    } else if(direction == "SW"){
        if(isOnBottomRow){
            if(isOnSWCorner){
                return processPerRow;
                //Process on the NE corner
            }
            if(isOnSECorner){
                return processPerRow - 1;
                //Process on the top row
            }
            return worldRank % processPerRow - 1;
        } else {
            if(isOnLeftmostColumn){
                return worldRank + 2 * processPerRow - 1;
                //Process on the top row
            }
            return worldRank + processPerRow - 1;
        }
    } else if(direction == "SE"){
        if(isOnBottomRow){
            if(isOnSECorner){
                return 1;
                //Process on the NW corner
            }
            return worldRank % processPerRow + 1;
        } else {
            if(isOnRightmostColumn){
                return worldRank + 1;
                //Process on the leftmost column
            }
            return worldRank + processPerRow + 1;
        }
    }
}

int* getPlace(string direction, int index){
    //Returns the pointer to the place of the
    //neighbor cell in given direction
    if(direction == "N"){
        return &ups[index];
    } else if(direction == "S") {
        return &downs[index];
    } else if(direction == "W"){
        return &lefts[index];
    } else if(direction == "E"){
        return &rights[index];
    } else if(direction == "NW"){
        return &corners[0];
    } else if(direction == "NE"){
        return &corners[1];
    } else if(direction == "SW"){
        return &corners[2];
    } else {
        return &corners[3];
    }
}

int getFromArrayOrReceive(string direction, int index){
    //Gets the neighbor of the cell, if it is not received
    //already, receives it and stores it in relevant array
    int result;//Neighbor cell gotten from neighbor process
    int *place = getPlace(direction, index);
    //Place of the desired cell in local arrays
    if(direction.length() == 2){
        index = 0;
        //Tag will be zero if the direction is diagonal
    }
    if(*place == NOT_RECEIVED){//Neighbor cell is not received yet
        MPI_Recv(&result, 1, MPI_INT, getWorldRank(direction), index, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //Receive the cell
        *place = result;
        //Store it in local array
    } else {
        result = *place;
        //Set the result from local array
    }
    return result;
}

int nOfLivingNeighbors(int i, int j){
    //Returns the number of living neighbors of the cell
    vector<vector<int>>* gameMap = stage % 2 == 1 ? &gameMaps[1] : &gameMaps[0];
    //Map to be used is chosen by looking at the stage
    int lastIndex = edgeSize - 1;//Last column index
    int N, NE, NW, W, E, S, SE, SW;//Neighbors in directions
    if(i == 0){//i is in the first row
        N = getFromArrayOrReceive("N", j);
        if(j == 0) {//j is in the first column
            NW = getFromArrayOrReceive("NW", j - 1);
        } else {
            NW = getFromArrayOrReceive("N", j - 1);
        }
        if(j == lastIndex) {//j is in the last column
            NE = getFromArrayOrReceive("NE", j + 1);
        } else {
            NE = getFromArrayOrReceive("N", j + 1);
        }
    } else {
        N = (*gameMap)[i - 1][j];
        if(j == 0){//j is in the first column
            NW = getFromArrayOrReceive("W", i - 1);
        } else {
            NW = (*gameMap)[i - 1][j - 1];
        }
        if(j == lastIndex){//j is in the last column
            NE = getFromArrayOrReceive("E", i - 1);
        } else {
            NE = (*gameMap)[i - 1][j + 1];
        }
    }
    if(i == lastIndex){//i is in the last row
        S = getFromArrayOrReceive("S", j);
        if(j == 0) {//j is in the first column
            SW = getFromArrayOrReceive("SW", j - 1);
        } else {
            SW = getFromArrayOrReceive("S", j - 1);
        }
        if(j == lastIndex) {//j is in the last column
            SE = getFromArrayOrReceive("SE", j + 1);
        } else {
            SE = getFromArrayOrReceive("S", j + 1);
        }
    } else {
        S = (*gameMap)[i + 1][j];
        if(j == 0){//j is in the first column
            SW = getFromArrayOrReceive("W", i + 1);
        } else {
            SW = (*gameMap)[i + 1][j - 1];
        }
        if(j == lastIndex){//j is in the last column
            SE = getFromArrayOrReceive("E", i + 1);
        } else {
            SE = (*gameMap)[i + 1][j + 1];
        }
    }
    if(j == lastIndex){//j is in the last column
        E = getFromArrayOrReceive("E", i);
    } else {
        E = (*gameMap)[i][j + 1];
    }
    if(j == 0){//j is in the first column
        W = getFromArrayOrReceive("W", i);
    } else {
        W = (*gameMap)[i][j - 1];
    }
    return N + NE + NW + S + SE + SW + E + W;
}

void update(int i, int j){
    int livingNeighboorCount = nOfLivingNeighbors(i, j);
    //Select the old game map and the new game map by looking at the stage
    vector<vector<int>>* newGameMap = stage % 2 == 1 ? &gameMaps[0] : &gameMaps[1];
    vector<vector<int>>* oldGameMap = stage % 2 == 1 ? &gameMaps[1] : &gameMaps[0];
    if(livingNeighboorCount < 2 || livingNeighboorCount > 3){
        (*newGameMap)[i][j] = 0;
        //Overpopulation or loneliness
    } else if(livingNeighboorCount == 3){
        (*newGameMap)[i][j] = 1;
        //Reproduce
    } else {
        (*newGameMap)[i][j] = (*oldGameMap)[i][j];
    }
}

void readInput(string filename){
    //Read input to the game map array
    if(worldRank == MASTER_PROCESS_RANK){
        gameMaps.push_back(vector<vector<int>>());
        string line;
        ifstream inputFile;
        inputFile.open(filename);
        while (getline(inputFile, line)){
            gameMaps[0].push_back(vector<int>());
            for(int i = 0; i < line.length() - 1; i+=2){
                gameMaps[0][gameMaps[0].size() - 1].push_back(line[i] - '0');
                //Store chars as int in gameMap array
            }
        }
        inputFile.close();
    }
}

void writeOutput(string filename){
    //Write the game map to the output file
    if(worldRank == MASTER_PROCESS_RANK){
        ofstream outputFile;
        outputFile.open(filename);
        for(int i = 0; i < gameMaps[0].size(); i++){
            for(int j = 0; j < gameMaps[0][i].size(); j++){
                outputFile << gameMaps[0][i][j] << " ";
            }
            outputFile << endl;
        }
        outputFile.close();
    }
}

void updateAll() {
    if(worldRank != MASTER_PROCESS_RANK){
        for(int i = 0; i < edgeSize; i++){
            for(int j = 0; j < edgeSize; j++){
                update(i, j);//Update all cells of game map
            }
        }
        //Reinitialize neighbor vectors
        ups = vector<int>(edgeSize, NOT_RECEIVED);
        downs = vector<int>(edgeSize, NOT_RECEIVED);
        lefts = vector<int>(edgeSize, NOT_RECEIVED);
        rights = vector<int>(edgeSize, NOT_RECEIVED);
        corners = vector<int>(4, NOT_RECEIVED);
        stage++;//Increment the stage
    }
}

void divideGame(){
    //Divide the game to the processes
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);//get world size
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);//get world rank
    if(worldRank == MASTER_PROCESS_RANK){
        edgeSize = gameMaps[0].size() / (int)sqrt(worldSize - 1);//set edge size
        for(int i = 1; i < worldSize; i++){
            int tag = 1;//Tag for cells
            int j = (i - 1) / (int)sqrt(worldSize - 1) * edgeSize;
            //Beginning of the row index
            MPI_Send(&edgeSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            //Send edge size to the process
            int boundaryJ = j + edgeSize;
            //Boundary of the row
            for(; j < boundaryJ; j++){
                int k = (i - 1) % (int)sqrt(worldSize - 1) * edgeSize;
                //Beginning of the column index
                int boundaryK = k + edgeSize;
                //Boundary of the column
                for(; k < boundaryK; k++){
                    //Send the cell and increment the tag
                    MPI_Send(&gameMaps[0][j][k], 1, MPI_INT, i, tag, MPI_COMM_WORLD);
                    tag++;
                }
            }
        }
    } else {
        //Initialize the maps
        gameMaps = vector<vector<vector<int>>>(2);
        //Receive the edge size
        MPI_Recv(&edgeSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //Set the tag for the cells
        int tag = 1;
        //Initalize the neighbor vectors
        ups = vector<int>(edgeSize, NOT_RECEIVED);
        downs = vector<int>(edgeSize, NOT_RECEIVED);
        lefts = vector<int>(edgeSize, NOT_RECEIVED);
        rights = vector<int>(edgeSize, NOT_RECEIVED);
        corners = vector<int>(4, NOT_RECEIVED);
        for(int i = 0; i < edgeSize; i++){
            //Push empty vectors to fill as rows
            gameMaps[0].push_back(vector<int>());
            gameMaps[1].push_back(vector<int>());
            for(int j=0; j < edgeSize; j++){
                //get new value of the cell
                int newVal;
                MPI_Recv(&newVal, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //Set the value as the cell of the arrays
                gameMaps[0][i].push_back(newVal);
                gameMaps[1][i].push_back(newVal);
                //Increment the tag for next cell
                tag++;
            }
        }
    }
}

void sendAll(){
    //Processes send the information needed to the adjacent cells
    if(worldRank != MASTER_PROCESS_RANK){
        //Get the current game map
        vector<vector<int>>* gameMap = stage % 2 == 1 ? &gameMaps[1] : &gameMaps[0];
        //Send the left and right columns
        for(int j = 0; j < edgeSize; j++){
            MPI_Send(&(*gameMap)[0][j], 1, MPI_INT, getWorldRank("N"), j, MPI_COMM_WORLD);
            MPI_Send(&(*gameMap)[edgeSize - 1][j], 1, MPI_INT, getWorldRank("S"), j, MPI_COMM_WORLD);
        }
        //Send the top and bottom rows
        for(int j = 0; j < edgeSize; j++){
            MPI_Send(&(*gameMap)[j][0], 1, MPI_INT, getWorldRank("W"), j, MPI_COMM_WORLD);
            MPI_Send(&(*gameMap)[j][edgeSize - 1], 1, MPI_INT, getWorldRank("E"), j, MPI_COMM_WORLD);
        }
        //Send the corners
        MPI_Send(&(*gameMap)[0][edgeSize - 1], 1, MPI_INT, getWorldRank("NE"), 0, MPI_COMM_WORLD);
        MPI_Send(&(*gameMap)[0][0], 1, MPI_INT, getWorldRank("NW"), 0, MPI_COMM_WORLD);
        MPI_Send(&(*gameMap)[edgeSize - 1][edgeSize - 1], 1, MPI_INT, getWorldRank("SE"), 0, MPI_COMM_WORLD);
        MPI_Send(&(*gameMap)[edgeSize - 1][0], 1, MPI_INT, getWorldRank("SW"), 0, MPI_COMM_WORLD);
    }
}

void recollect(){
    if(worldRank != MASTER_PROCESS_RANK){
        //Get the current game map
        vector<vector<int>> gameMap = *(stage % 2 == 1 ? &gameMaps[1] : &gameMaps[0]);
        int tag = 0;//Tag for the cell
        for(int i = 0; i < edgeSize; i++){
            for(int j = 0; j < edgeSize; j++){
                //Send each cell
                MPI_Send(&gameMap[i][j], 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
                tag++;
            }
        }
    } else {
        for(int p = 1; p < worldSize; p++) {
            //Tag for the cells
            int tag = 0;
            for(int i = 0; i < edgeSize; i++) {
                for(int j = 0; j < edgeSize; j++) {
                    //Get the indexes of the cell
                    int indexI = (p - 1) / (int)sqrt(worldSize - 1) * edgeSize + i;
                    int indexJ = (p - 1) % (int)sqrt(worldSize - 1) * edgeSize + j;
                    //Receive the cell data
                    MPI_Recv(&gameMaps[0][indexI][indexJ], edgeSize, MPI_INT, p, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);    
                    tag++;
                }
            }
        }
    }
}

int main(int argc, char *args[]){
    string inputFileName = args[1];
    string outputFileName = args[2];
    int nOfIterations = stoi(args[3]);
    MPI_Init(NULL, NULL);
    readInput(inputFileName);
    divideGame();
    for(int i = 0; i < nOfIterations; i++){
        sendAll();
        updateAll();        
    }
    recollect();
    writeOutput(outputFileName);
    MPI_Finalize();
}
