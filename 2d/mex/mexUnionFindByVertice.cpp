/*==========================================================
 * mexUnionFindBodies.cpp 
 * To compile type: mex -R2018a mexUnionFindByVertice.cpp 
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <math.h>
#include <unordered_map>

class UF { 
  public:
  int *id, setCount, *sizeOfSet, *claimedVerts;
// Create an empty union find data structure with N isolated sets.
UF(int numElements, int numVertice) {
    setCount = 0; 
    id = new int[numElements]; 
    sizeOfSet = new int[numElements]; 
    claimedVerts = new int[numVertice];

    //initially things point towards -1 until they are inserted
    for (int i = 0; i<numElements; ++i) {
        id[i] = -1; 
        sizeOfSet[i] = 1; 
    }

    for (int i = 0; i<numVertice; ++i) { 
        claimedVerts[i] = -1; 
    }
}
~UF() { delete[] id; delete[] sizeOfSet; delete[] claimedVerts;}

//enables the creation of sets for element
bool insert(size_t p, bool * forcedElasticElements){
    if (!forcedElasticElements[p])
    {
        id[p] = p;
        setCount += 1;
        return true;
    }
    return false;
}

// Return the id of component corresponding to object p.
int find(int p) {
    if (p == -1 || id[p] == -1)
        return -1;

    int root = p;
    //find root
    while (root != id[root])    
        root = id[root];
    
    //compress all the path to root
    while (p != root) { 
        int newp = id[p]; 
        id[p] = root; 
        p = newp; 
    }
    return root;
}
// Replace sets containing x and y with their union.
void merge(int x, int y) {
    int i = find(x); 
    int j = find(y); 

    if (i == j || i == -1 || j == -1) {return;}
        
    // make smaller root point to larger one
    if (sizeOfSet[i] < sizeOfSet[j]) { 
        id[i] = j;
        sizeOfSet[j] += sizeOfSet[i]; 
    }
    else { 
        id[j] = i;
        sizeOfSet[i] += sizeOfSet[j]; 
    }
    setCount--;
}
// Are objects x and y in the same set?
bool connected(int x, int y) { return find(x) == find(y); }
// Return the number of disjoint sets.
int count() { return setCount; }
};

void copyLayer(double *rigidBodySetsPerLayer, const int numElements, size_t layerIndex, UF* rigidSets, double* numRigids){
    int counter = 1;//starts at 1 here because matlab is 1-indexed
    //map element roots to rigid body ID
    std::unordered_map<int, int> mapToRigid;
    
    numRigids[layerIndex] = rigidSets->count();
    for (size_t i = 0; i < numElements; ++i){
        //check if the element is elastic or rigid
        size_t position = i + numElements*layerIndex;
        if (rigidSets->id[i] == -1){
            rigidBodySetsPerLayer[position] = 0;
        } else {
            //map from find to counter
            int root = rigidSets->find(i);
            //if exists, then element has rigid id = counter
            auto rigidID = mapToRigid.find(root);
            if (rigidID != mapToRigid.end()){ 
                rigidBodySetsPerLayer[position] = rigidID->second;
            } else{
                //if not then add entry to map(find, counter) and increase counter
                mapToRigid.insert({root,counter});
                rigidBodySetsPerLayer[position] = counter;
                counter += 1;
            }
        }
    }
}

/* 
 * The gateway function
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    //assumes the layers are in ascending order and not repeating (please filter before sending data)
    //assumes that the first layer passed is always 0
    double *elementsInStrainRateOrder = mxGetDoubles(prhs[0]);
    size_t numElements = mxGetM(prhs[0]);
    double *layers = mxGetDoubles(prhs[1]);
    const int numLayers= mxGetN(prhs[1]);
    bool* forcedElasticElements = mxGetLogicals(prhs[2]);
    int* T   =(int*)mxGetData(prhs[3]); 
    const int vertsPerT= mxGetN(prhs[3]);
    double numVerts   = *mxGetDoubles(prhs[4]);

    bool isLogicalforcedElasticElements = mxIsLogical(prhs[2]);
    if (!isLogicalforcedElasticElements) {
        mexErrMsgIdAndTxt("ARP:mexUnionFindBodies:type", "forcedElasticElements must be logical");
    }
    if (numElements != mxGetM(prhs[2])) {
        mexErrMsgIdAndTxt("ARP:mexUnionFindBodies:size", "forcedElasticElements: Expecting a column vector of numElements");
    }

    if ( nrhs != 5 ) {
        mexErrMsgIdAndTxt("ARP:mexUnionFindBodies:nrhs","5 inputs required.");
    }
    if ( nlhs != 2 ) {
        mexErrMsgIdAndTxt("ARP:mexUnionFindBodies:nlhs","2 outputs required.");
    }

    if( mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("ARP:mexUnionFindBodies:triVertIDint32","triVertIDint32 must be a int matrix.");
    }

    plhs[0] = mxCreateDoubleMatrix( numElements, numLayers, mxREAL );
    double *rigidBodySetsPerLayer = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix( numLayers, 1, mxREAL );
    double *numRigids = mxGetDoubles(plhs[1]);
    
    //initialize the union-find
    UF* rigidSets = new UF(numElements, (int)floor(numVerts));
    copyLayer(rigidBodySetsPerLayer, numElements, 0, rigidSets, numRigids);

    //for all layers
    int layerStart = 0;
    for(size_t layer = 0; layer<numLayers; ++layer){
        int layerEnd = floor((layers[layer]/100.0)*numElements);
        //insert elements of elementsInStrainRateOrder from lastLayerEnd to layerEnd
        for(size_t e = layerStart; e < layerEnd; ++e){
            int element = elementsInStrainRateOrder[e]-1; //-1 here because matlab is 1-indexed
            bool successfulInsertion = rigidSets->insert(element,forcedElasticElements);

            if(successfulInsertion){
                for(size_t vertOfT = 0; vertOfT<vertsPerT; ++vertOfT){
                    int vert = T[element + numElements*vertOfT]-1;
                    rigidSets->merge(element,rigidSets->claimedVerts[vert]);
                    rigidSets->claimedVerts[vert] = element;
                }
            }
        }
        copyLayer(rigidBodySetsPerLayer,numElements,layer,rigidSets, numRigids);
        layerStart = layerEnd;
    }
    delete rigidSets;
}
