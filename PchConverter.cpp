/* -------------------------------*\
| Copyright (C) 2020 Pavel Lakiza  |
\* -------------------------------*/

// The command to perform export from Femap: PARAM,EXTOUT,DMIGPCH

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <unordered_map>
#include <iomanip>    
#include <algorithm>
#include <numeric>

using std::vector;
using std::string;
using std::cout;
using std::cin;
using std::endl;

using IndexType = size_t;

// A hash function used to hash a pair of any kind 
struct hash_pair {
	template <class T1, class T2>
	size_t operator()(const std::pair<T1, T2>& p) const
	{
		auto hash1 = std::hash<T1>{}(p.first);
		auto hash2 = std::hash<T2>{}(p.second);
		return hash1 ^ hash2;
	}
};

template <typename T>
void releaseVector(vector<T>& vec) {
	vec.clear();
	vec.shrink_to_fit();
}

enum class MatrixType { Stiffness, Mass, Damping, Mapping, Null = -1 };

// A base matrix class
class StructuralMatrix {
public:
	StructuralMatrix(string const& name, string const& outputFileName);
	~StructuralMatrix();
	// Info methods
	bool isEmpty() const { return realSize_ == 0; } // Check whether the matrix empty or not
	IndexType size() const { return realSize_; } // Get the size of the matrix
	string const& name() const { return name_; } // Get the name of the matrix
	// Managing methods
	void specifySize(IndexType size); // Remember size for further allocation
	void resize(); // Change the size of the matrix
	void sortByColumn(); // Sort the matrix by the column index
	// In-out methods
	int write(unsigned short outPrecision) const; // Write the matrix
	// Get the elements of the matrix
	IndexType& firstInd(IndexType index) { return firstInd_[index]; }
	IndexType& secondInd(IndexType index) { return secondInd_[index]; }
	double& value(IndexType index) { return values_[index]; }
private:
	string name_;
	string outputFileName_ = "mat.prn";
	IndexType realSize_ = 0;
	IndexType specifiedSize_ = 0;
	vector<IndexType> firstInd_;
	vector<IndexType> secondInd_;
	vector<double> values_;
};

// Constructor
StructuralMatrix::StructuralMatrix(string const& name, string const& outputFileName) :
	name_(name), 
	outputFileName_(outputFileName)
{

}

// Destructor

StructuralMatrix::~StructuralMatrix() {
	realSize_ = 0;
	specifiedSize_ = 0;
	releaseVector(firstInd_);
	releaseVector(secondInd_);
	releaseVector(values_);
}

// Remember size for further allocation
void StructuralMatrix::specifySize(IndexType size) {
	specifiedSize_ = size;
}

// Changing the size of the structural matrix
void StructuralMatrix::resize() {
	if (specifiedSize_ <= 0)
		return;
	realSize_ = specifiedSize_;
	firstInd_.resize(realSize_);
	secondInd_.resize(realSize_);
	values_.resize(realSize_);
}

// Sorting the structural matrix by the column index
void StructuralMatrix::sortByColumn() {
	vector<IndexType> sortIndexes(realSize_);
	std::iota(sortIndexes.begin(), sortIndexes.end(), 0);
	// Copying the current vectors in order to sort them later
	vector<IndexType> vecFirst = firstInd_;
	vector<IndexType> vecSecond = secondInd_;
	vector<double> vecValues = values_;
	// Sorting by the column index 
	std::stable_sort(sortIndexes.begin(), sortIndexes.end(), [&vecSecond](IndexType i, IndexType j) { return vecSecond[i] < vecSecond[j]; });
	IndexType tInd = 0;
	for (IndexType i = 0; i != realSize_; ++i) {
		tInd = sortIndexes[i];
		firstInd_[i] = vecFirst[tInd];
		secondInd_[i] = vecSecond[tInd];
		values_[i] = vecValues[tInd];
	}
}

// Get the type of the strutural matrix
MatrixType resolveType(string const& strType) {
	if ( !strType.compare("KAAX") ) return MatrixType::Stiffness;
	if ( !strType.compare("MAAX") ) return MatrixType::Mass;
	if ( !strType.compare("DAAX") ) return MatrixType::Damping;
	if ( !strType.compare("VAX") ) return MatrixType::Mapping;
	return MatrixType::Null;
}

// Open file for reading
bool openFileWithMessages(string& fileName, std::ifstream& stream) {
	stream.open(fileName, std::ios::in);
	bool isOpened = stream.is_open();
	if (isOpened)
		cout << "* '" << fileName << "' was successfully opened" << endl;
	else
		cout << "* An error occured while opening '" << fileName << "'" << endl;
	return isOpened;
}

// Organize a text dialog with a user to open a file
void organizeOpenFileDialog(string& fileName, std::ifstream& filePch, string const& message) {
	while (true) {
		cout << message;
		cin >> fileName;
		openFileWithMessages(fileName, filePch);
	}
}

// Organizing a text dialog to get a number
template < typename T >
T organizeNumberDialog(string const& message, T const& hintValue) {
	T value;
	while ( true ) {
		cout << message;
		cin >> value;
		if ( !cin.fail() )
			break;
	}
	return value;
}

// Writing the structural matrix
int StructuralMatrix::write(unsigned short outPrecision) const {
	std::ofstream file;
	file.open(outputFileName_, std::ios::out);
	if (!file.is_open()) {
		cout << "An error occurred while saving '" << outputFileName_ << "'" << endl;
		return -1;
	}
	file.precision(outPrecision);
	for (IndexType i = 0; i != realSize_; ++i) {
		file << std::fixed << std::setw(outPrecision) << firstInd_[i] << std::setw(outPrecision) << secondInd_[i] << ' ';
		file << std::scientific << '\t' << values_[i] << endl;
	}
	file.close();
	return 0;
}

// Writing the mapping matrix
int writeMappingMatrix(string const& fileName, IndexType size, vector<IndexType> const& vecEquNodes, vector<IndexType> const& vecEquDofs, unsigned short outPrecision) {
	std::ofstream file;
	file.open(fileName, std::ios::out);
	if ( !file.is_open() ) {
		cout << "An error occurred while saving '" << fileName << "'" << endl;
		return -1;
	}
	vector<string> const dofsNames = { "UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ" };
	file << std::setw(outPrecision) << "Matrix Eqn" << std::setw(outPrecision) << "Node" << std::setw(outPrecision) << "DOF" << endl;
	for (IndexType i = 0; i != size; ++i) {
		file << std::setw(outPrecision) << i + 1;
		file << std::setw(outPrecision) << vecEquNodes[i] << std::setw(outPrecision) << dofsNames[vecEquDofs[i] - 1] << endl;
	}
	file.close();
	return 0;
}

// Writing the nodes file
int writeNodes(string const& fileName, vector<IndexType> const& indexes, vector<double> const& XCoord, vector<double> const& YCoord, vector<double> const& ZCoord, unsigned short outPrecision) {
	std::ofstream file;
	file.open(fileName, std::ios::out);
	if (!file.is_open()) {
		cout << "An error occurred while saving '" << fileName << "'" << endl;
		return -1;
	}
	IndexType size = indexes.size();
	for (IndexType i = 0; i != size; ++i) {
		file << std::fixed << std::setw(outPrecision) << indexes[i];
		file << std::scientific << '\t' << XCoord[i] << '\t' << YCoord[i] << '\t' << ZCoord[i] << endl;
	}
	file.close();
	return 0;
}

int main(){

	// Defining reading constants
	const unsigned short LENGTH_MATRIX_HEADER = 6;
	const unsigned short NUMBER_OF_SKIPPED_SYMBOLS = 16384;
	const unsigned short OUTPUT_PRECISION = 10;
	const string FILE_LOADER_NAME = "Loader.txt";
		
	string fileName;
	std::ifstream filePch, fileDat;
	double geometryScaleFactor = 1.0;
	bool isRead = false;

	// Reading input files according to a loader
	std::ifstream fileLoader(FILE_LOADER_NAME, std::ios::in);
	isRead = fileLoader.is_open();
	if (isRead) {
		cout << "Reading files according to the file loader:" << endl;
		// Punch file
		fileLoader >> fileName;
		isRead = isRead && openFileWithMessages(fileName, filePch);
		// Analysis input file
		fileLoader >> fileName;
		isRead = isRead && openFileWithMessages(fileName, fileDat);
		// Scale
		fileLoader >> geometryScaleFactor;
		if (!isRead)
			cout << "The format of the file loader is not correct. Switching to user input..." << endl;
	} 

	// Reading input files from the keyboard
	if (!isRead){
		// Opening an input Nastran punch file
		organizeOpenFileDialog(fileName, filePch, "Specify the path of the Nastran punch file (with '.pch'): ");
		// Opening an analysis input file
		organizeOpenFileDialog(fileName, fileDat, "Specify the path of the analysis input file (with '.dat'): ");
		geometryScaleFactor = organizeNumberDialog("Specify the geometry scale factor: ", 1.0);
	}

	auto startTime = std::chrono::steady_clock::now();

	// Extracting the nodal coordinates from the analysis input file
	cout << "Extracting the nodal coordinates from the analysis input file...";
	string tempString;
		// Counting the nodes number
	IndexType nNodes = 0;
	bool isGrid = false;
	while ( !fileDat.eof() ) {
		fileDat >> tempString;
		if (tempString.find("GRID") != std::string::npos) {
			++nNodes;
			isGrid = true;
			if (tempString[tempString.size() - 1] == '*')
				fileDat.ignore(NUMBER_OF_SKIPPED_SYMBOLS, '\n');
		} else if (isGrid) {
			break;
		}
		fileDat.ignore(NUMBER_OF_SKIPPED_SYMBOLS, '\n');
	}
		// Reading the nodes coordinates
	vector<IndexType> nodesNumbers(nNodes);
	vector<double> nodesX(nNodes);
	vector<double> nodesY(nNodes);
	vector<double> nodesZ(nNodes);
	nNodes = 0;
	isGrid = false;
	fileDat.clear();
	fileDat.seekg(0);
	IndexType j;
	while (!fileDat.eof()) {
		fileDat >> tempString;
		if (tempString.find("GRID") != std::string::npos) {
			fileDat >> nodesNumbers[nNodes] >> j >> nodesX[nNodes] >> nodesY[nNodes];
			if (tempString[tempString.size() - 1] == '*') {
				fileDat.ignore(NUMBER_OF_SKIPPED_SYMBOLS, '\n');
				fileDat >> tempString;
			}
			fileDat >> nodesZ[nNodes];
			// Scaling 
			nodesX[nNodes] *= geometryScaleFactor;
			nodesY[nNodes] *= geometryScaleFactor;
			nodesZ[nNodes] *= geometryScaleFactor;
			++nNodes;
			isGrid = true;
		} else if (isGrid) {
			break;
		}
		fileDat.ignore(NUMBER_OF_SKIPPED_SYMBOLS, '\n');
	}
	fileDat.close();
	cout << "OK" << endl;

	// Writing the nodes file
	cout << "Writing the nodes file...";
	IndexType i;
	i = writeNodes("nodes2.prn", nodesNumbers, nodesX, nodesY, nodesZ, OUTPUT_PRECISION);
	if (i == 0) cout << "OK" << endl;
	// Releasing nodes
	releaseVector(nodesNumbers);
	releaseVector(nodesX);
	releaseVector(nodesY);
	releaseVector(nodesZ);

	// Announcing of the structural matrices
	StructuralMatrix stiffnessMatrix("stiffness", "matK.prn");
	StructuralMatrix massMatrix("mass", "matM.prn");
	StructuralMatrix dampingMatrix("damping", "matD.prn");

	// Counting the numbers of nonzero elements of the matrices
	cout << "Counting the numbers of nonzero elements of the matrices...";
	bool isMatrixChanged = false;
	StructuralMatrix* ptrMatrix = nullptr;
	IndexType nNonzeroElements = 0;
	IndexType nMatrices = 0;
	bool isExit = false;
	while (!filePch.eof()) {
		filePch >> tempString; // DMIG or DMIG* or *
		// Check if the matrix was changed or not
		isMatrixChanged = false;
		if (!tempString.compare("DMIG")) {
			// Reading the type of the matrix
			isMatrixChanged = true;
			filePch >> tempString;
			// Resizing the matrices
			if (ptrMatrix != nullptr)
				ptrMatrix->specifySize(nNonzeroElements); // (!) Memory is not allocated here
			// Check the type of the matrix
			switch ( resolveType(tempString) ) {
			case MatrixType::Stiffness:
				ptrMatrix = &stiffnessMatrix;
				break;
			case MatrixType::Mass:
				ptrMatrix = &massMatrix;
				break;
			case MatrixType::Damping:
				ptrMatrix = &dampingMatrix;
				break;
			case MatrixType::Mapping:
				isExit = true;
				break;
			case MatrixType::Null:
				ptrMatrix = nullptr;
				continue;
			}
			if (isExit)
				break;
			// Getting the size of matrices
			for (unsigned short i = 0; i != LENGTH_MATRIX_HEADER - 2; ++i)
				filePch >> tempString;
			filePch >> nMatrices; // Size of the structural matrices (free model)
			nNonzeroElements = 0;
		} else if (ptrMatrix != nullptr && !tempString.compare("*")){
			++nNonzeroElements;
		}
	}
	cout << "OK" << endl;

	// Reading the mapping matrix
	cout << "Reading the mapping matrix...";
	vector<IndexType> vecEquNodes(nMatrices);
	vector<IndexType> vecEquDofs(nMatrices);
	std::unordered_map<std::pair<IndexType, IndexType>, IndexType, hash_pair> mapEqu; // Associative container which is used to get an equation number by a node number and its dof  
	for (i = 0; i != 2; ++i)
		filePch.ignore(NUMBER_OF_SKIPPED_SYMBOLS, '\n');
	nMatrices = 0; // Real number of dofs
	double tempValue;
	while ( true ) {
		filePch >> tempString;
		if (filePch.eof() || !tempString.compare("DMIG"))
			break;
		filePch >> vecEquNodes[nMatrices] >> vecEquDofs[nMatrices] >> tempValue;
		mapEqu.emplace(std::make_pair(vecEquNodes[nMatrices], vecEquDofs[nMatrices]), nMatrices + 1);
		++nMatrices;
	}
	cout << "OK" << endl;

	// Clearance of working variables
	ptrMatrix = nullptr;
	isExit = false;

	// Filling the structural matrices
	filePch.clear();
	filePch.seekg(0);
	IndexType iEqu = 0, jEqu = 0;
	while (!filePch.eof()) {
		filePch >> tempString; // DMIG or DMIG* or *
		// Check if the matrix was changed or not
		isMatrixChanged = false;
		if (!tempString.compare("DMIG")) {
			// Reading the type of the matrix
			isMatrixChanged = true;
			filePch >> tempString;
			// Sorting and writing the previous matrix
			if (ptrMatrix != nullptr) {
				cout << "OK" << endl;
				string const& matrixName = ptrMatrix->name();
				cout << "Sorting the " << matrixName << " matrix by column index...";
				ptrMatrix->sortByColumn();
				cout << "OK" << endl;
				cout << "Writing the " << matrixName << " matrix...";
				i = ptrMatrix->write(OUTPUT_PRECISION);
				if (i == 0) cout << "OK" << endl;
				ptrMatrix->~StructuralMatrix();
			}
			// Check the type of the matrix
			switch ( resolveType(tempString) ) {
			case MatrixType::Stiffness:
				ptrMatrix = &stiffnessMatrix;
				break;
			case MatrixType::Mass:
				ptrMatrix = &massMatrix;
				break;
			case MatrixType::Damping:
				ptrMatrix = &dampingMatrix;
				break;
			case MatrixType::Mapping:
				isExit = true;
				break;
			case MatrixType::Null:
				ptrMatrix = nullptr;
				continue;
			}
			if (isExit)
				break;
			ptrMatrix->resize();
			cout << "Filling the " << ptrMatrix->name() << " matrix...";
			filePch.ignore(NUMBER_OF_SKIPPED_SYMBOLS, '\n');
			nNonzeroElements = 0; 
		} else if (ptrMatrix != nullptr) {
			if (!tempString.compare("DMIG*")) {	   // Column Index
				filePch >> tempString >> i >> j;
				jEqu = mapEqu[std::make_pair(i, j)];
			} else if (!tempString.compare("*")) { // Row index
				filePch >> i >> j >> tempValue;
				iEqu = mapEqu[std::make_pair(i, j)];
				// Filling the current matrix (swap indexes: lower triangle -> upper triangle)
				ptrMatrix->firstInd(nNonzeroElements) = jEqu;
				ptrMatrix->secondInd(nNonzeroElements) = iEqu;
				ptrMatrix->value(nNonzeroElements) = tempValue;
				++nNonzeroElements;
			}
		}
	}
	filePch.close();

	// Writing the mapping and structural matrices
	cout << "Writing the mapping matrix...";
	i = writeMappingMatrix("matK.mapping", nMatrices, vecEquNodes, vecEquDofs, OUTPUT_PRECISION);
	if ( i == 0 ) cout << "OK" << endl;
		
	// Evaluation of the duration of the conversion
	auto endTime = std::chrono::steady_clock::now();
	std::chrono::duration<double> duration = endTime - startTime;
	cout << "Duration: " << duration.count() << " s\n";

	// Pause
	cout << "Press any key to continue . . .";
	cin.ignore();
	cin.get();
	return 0;
}