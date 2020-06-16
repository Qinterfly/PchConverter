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

enum class MatrixType { Stiffness, Mass, Damping, Mapping, Null = -1 };

// A base matrix class
class StructuralMatrix {
public:
	StructuralMatrix() = default;
	~StructuralMatrix() = default;
	// Info methods
	bool isEmpty() const { return size_ == 0; } // Check whether the matrix empty or not
	int size() const { return size_; } // Get the size of the matrix
	// Managing methods
	void resize(long size); // Change the size of the matrix
	void sortByColumn(); // Sort the matrix by the column index
	// In-out methods
	int write(string const& fileName, unsigned short outPrecision) const; // Write the matrix
	// Get the elements of the matrix
	int& firstInd(int index) { return firstInd_[index]; }
	int& secondInd(int index) { return secondInd_[index]; }
	double& value(int index) { return values_[index]; }
private:
	int size_ = 0;
	vector<int> firstInd_;
	vector<int> secondInd_;
	vector<double> values_;
};

// Changing the size of the structural matrix
void StructuralMatrix::resize(long size) {
	size_ = size;
	firstInd_.resize(size_);
	secondInd_.resize(size_);
	values_.resize(size_);
}

// Sorting the structural matrix by the column index
void StructuralMatrix::sortByColumn() {
	vector<int> sortIndexes(size_);
	std::iota(sortIndexes.begin(), sortIndexes.end(), 0);
	// Copying the current vectors in order to sort them later
	vector<int> vecFirst = firstInd_;
	vector<int> vecSecond = secondInd_;
	vector<double> vecValues = values_;
	// Sorting by the column index 
	std::stable_sort(sortIndexes.begin(), sortIndexes.end(), [&vecSecond](int i, int j) { return vecSecond[i] < vecSecond[j]; });
	int tInd = 0;
	for (int i = 0; i != size_; ++i) {
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

// Organizing a text dialog with a user to open a file
void organizeOpenFileDialog(string& fileName, std::ifstream& filePch, string const& message) {
	while (true) {
		cout << message;
		cin >> fileName;
		filePch.open(fileName, std::ios::in);
		if (filePch.is_open()) {
			cout << "* '" << fileName << "' was successfully opened" << endl;
			break;
		}
		else {
			cout << "* An error occured while opening '" << fileName << "'" << endl;
		}
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
int StructuralMatrix::write(string const& fileName, unsigned short outPrecision) const {
	std::ofstream file;
	file.open(fileName, std::ios::out);
	if (!file.is_open()) {
		cout << "An error occurred while saving '" << fileName << "'" << endl;
		return -1;
	}
	file.precision(outPrecision);
	for (int i = 0; i != size_; ++i) {
		file << std::fixed << std::setw(outPrecision) << firstInd_[i] << std::setw(outPrecision) << secondInd_[i] << ' ';
		file << std::scientific << '\t' << values_[i] << endl;
	}
	file.close();
	return 0;
}

// Writing the mapping matrix
int writeMappingMatrix(string const& fileName, int size, vector<int> const& vecEquNodes, vector<int> const& vecEquDofs, unsigned short outPrecision) {
	std::ofstream file;
	file.open(fileName, std::ios::out);
	if ( !file.is_open() ) {
		cout << "An error occurred while saving '" << fileName << "'" << endl;
		return -1;
	}
	vector<string> const dofsNames = { "UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ" };
	file << std::setw(outPrecision) << "Matrix Eqn" << std::setw(outPrecision) << "Node" << std::setw(outPrecision) << "DOF" << endl;
	for (int i = 0; i != size; ++i) {
		file << std::setw(outPrecision) << i + 1;
		file << std::setw(outPrecision) << vecEquNodes[i] << std::setw(outPrecision) << dofsNames[vecEquDofs[i] - 1] << endl;
	}
	file.close();
	return 0;
}

// Writing the nodes file
int writeNodes(string const& fileName, vector<int> const& indexes, vector<double> const& XCoord, vector<double> const& YCoord, vector<double> const& ZCoord, unsigned short outPrecision) {
	std::ofstream file;
	file.open(fileName, std::ios::out);
	if (!file.is_open()) {
		cout << "An error occurred while saving '" << fileName << "'" << endl;
		return -1;
	}
	int size = indexes.size();
	for (int i = 0; i != size; ++i) {
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
 
	// Opening an input Nastran punch file
	string fileName;
	std::ifstream filePch;
	organizeOpenFileDialog(fileName, filePch, "Specify the path of the Nastran punch file (with '.pch'): ");
	
	// Opening an analysis input file
	std::ifstream fileDat;
	double geometryScaleFactor = 1.0;
	organizeOpenFileDialog(fileName, fileDat, "Specify the path of the analysis input file (with '.dat'): ");
	geometryScaleFactor = organizeNumberDialog("Specify the geometry scale factor: ", 1.0);

	auto startTime = std::chrono::steady_clock::now();

	// Extracting the nodal coordinates from the analysis input file
	cout << "Extracting the nodal coordinates from the analysis input file...";
	string tempString;
		// Counting the nodes number
	int nNodes = 0;
	bool isGrid = false;
	while ( !fileDat.eof() ) {
		fileDat >> tempString;
		if (!tempString.compare("GRID")) {
			++nNodes;
			isGrid = true;
		} else if (isGrid) {
			break;
		}
		fileDat.ignore(NUMBER_OF_SKIPPED_SYMBOLS, '\n');
	}
		// Reading the nodes coordinates
	vector<int> nodesNumbers(nNodes);
	vector<double> nodesX(nNodes);
	vector<double> nodesY(nNodes);
	vector<double> nodesZ(nNodes);
	nNodes = 0;
	isGrid = false;
	fileDat.clear();
	fileDat.seekg(0);
	int j;
	while (!fileDat.eof()) {
		fileDat >> tempString;
		if (!tempString.compare("GRID")) {
			fileDat >> nodesNumbers[nNodes] >> j >> nodesX[nNodes] >> nodesY[nNodes] >> nodesZ[nNodes];
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

	// Announcing of the structural matrices
	StructuralMatrix stiffnessMatrix;
	StructuralMatrix massMatrix;
	StructuralMatrix dampingMatrix;

	// Counting the numbers of nonzero elements of the matrices
	cout << "Counting the numbers of nonzero elements of the matrices...";
	bool isMatrixChanged = false;
	StructuralMatrix* ptrMatrix = nullptr;
	int nNonzeroElements = 0; 
	int nMatrices = 0;
	bool isExit = false;
	while ( !filePch.eof() ) {
		filePch >> tempString; // DMIG or DMIG* or *
		// Check if the matrix was changed or not
		isMatrixChanged = false;
		if (!tempString.compare("DMIG")) {
			// Reading the type of the matrix
			isMatrixChanged = true;
			filePch >> tempString;
			// Resizing the matrices
			if (ptrMatrix != nullptr)
				ptrMatrix->resize(nNonzeroElements);
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
				break;
			}
			if (isExit)
				break;
			// Getting the size of matrices
			for (int i = 0; i != LENGTH_MATRIX_HEADER - 2; ++i)
				filePch >> tempString;
			filePch >> nMatrices; // Size of the structural matrices (free model)
			nNonzeroElements = 0;
		} else if ( !tempString.compare("*") ){
			++nNonzeroElements;
		}
	}
	cout << "OK" << endl;

	// Reading the mappnig matrix
	cout << "Reading the mapping matrix...";
	vector<int> vecEquNodes(nMatrices);
	vector<int> vecEquDofs(nMatrices);
	std::unordered_map<std::pair<int, int>, int, hash_pair> mapEqu; // Associative container which is used to get the equation number by the node number and its dof  
	int i;
	for (i = 0; i != 2; ++i)
		filePch.ignore(NUMBER_OF_SKIPPED_SYMBOLS, '\n');
	nMatrices = 0; // Real number of dofs
	double tempValue;
	while ( true ) {
		filePch >> tempString;
		if ( filePch.eof() || !tempString.compare("DMIG") )
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
	cout << "Filling the structural matrices...";
	filePch.clear();
	filePch.seekg(0);
	int iEqu = 0, jEqu = 0;
	while ( !filePch.eof() ) {
		filePch >> tempString; // DMIG or DMIG* or *
		// Check if the matrix was changed or not
		isMatrixChanged = false;
		if ( !tempString.compare("DMIG") ) {
			// Reading the type of the matrix
			isMatrixChanged = true;
			filePch >> tempString;
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
				break;
			}
			if ( isExit )
				break;
			// Getting the size of matrices
			for (i = 0; i != LENGTH_MATRIX_HEADER - 2; ++i)
				filePch >> tempString;
			nNonzeroElements = 0; 
		} else if ( !tempString.compare("DMIG*") ) { // Column index
			filePch >> tempString >> i >> j;
			jEqu = mapEqu[std::make_pair(i, j)];
		} else if ( !tempString.compare("*") ) { // Row index
			filePch >> i >> j >> tempValue;
			iEqu = mapEqu[std::make_pair(i, j)];
			// Filling the current matrix (swap indexes: lower triangle -> upper triangle)
			ptrMatrix->firstInd(nNonzeroElements) = jEqu;
			ptrMatrix->secondInd(nNonzeroElements) = iEqu;
			ptrMatrix->value(nNonzeroElements) = tempValue;
			++nNonzeroElements;
		}
	}
	filePch.close();
	cout << "OK" << endl;

	// Sorting the structural matrices by the column index
	cout << "Sorting the matrices by column index...";
	stiffnessMatrix.sortByColumn(); // Stiffness
	massMatrix.sortByColumn(); // Mass
	if ( !dampingMatrix.isEmpty() )
		dampingMatrix.sortByColumn(); // Damping
	cout << "OK" << endl;

	// Writing the mapping and structural matrices
	cout << "Writing the mapping and the structural matrices...";
	i = stiffnessMatrix.write("matK.prn", OUTPUT_PRECISION);
	i = massMatrix.write("matM.prn", OUTPUT_PRECISION);
	if ( !dampingMatrix.isEmpty() )
		i = dampingMatrix.write("matD.prn", OUTPUT_PRECISION);
	i = writeMappingMatrix("matK.mapping", nMatrices, vecEquNodes, vecEquDofs, OUTPUT_PRECISION);
	if ( i == 0 ) cout << "OK" << endl;
		
	// Writing the nodes file
	cout << "Writing the nodes file...";
	i = writeNodes("nodes2.prn", nodesNumbers, nodesX, nodesY, nodesZ, OUTPUT_PRECISION);
	if (i == 0) cout << "OK" << endl;

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