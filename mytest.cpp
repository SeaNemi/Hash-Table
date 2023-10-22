// CMSC 341 - Spring 2023 - Project 4
#include "vdetect.h"
#include <random>
#include <vector>
#include <string>
enum RANDOM {UNIFORMINT, UNIFORMREAL, NORMAL};
class Random {
public:
    Random(int min, int max, RANDOM type=UNIFORMINT, int mean=50, int stdev=20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL){
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            //the mean and standard deviation can change by passing new values to constructor 
            m_normdist = std::normal_distribution<>(mean,stdev);
        }
        else if (type == UNIFORMINT) {
            //the case of UNIFORMINT to generate integer numbers
            // Using a fixed seed value generates always the same sequence
            // of pseudorandom numbers, e.g. reproducing scientific experiments
            // here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min,max);
        }
        else{ //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
        }
    }
    void setSeed(int seedNum){
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }

    int getRandNum(){
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if(m_type == NORMAL){
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while(result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT){
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum(){
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result*100.0)/100.0;
        return result;
    }
    
    private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;//normal distribution
    std::uniform_int_distribution<> m_unidist;//integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal;//real uniform distribution

};

//Tester class
class Tester{
    public:
        //testCondition displays the output of whether a test passed or failed
        void testCondition(bool var);

        bool testInsertion(VDetect& vdetect);
        bool collisionInsertion(VDetect& vdetect);

        bool errorGetVirus(VDetect& vdetect);
        bool normalGetVirus(VDetect& vdetect);
        bool collisionGetVirus(VDetect& vdetect);

        bool normalRemove(VDetect& vdetect);
        bool collisionRemove(VDetect& vdetect);

        bool insertRehash(VDetect& vdetect);
        bool completeLambdaRehash(VDetect& vdetect);

        bool removeRehash(VDetect& vdetect);
        bool completeDeleteRehash(VDetect& vdetect);
};

unsigned int hashCode(const string str);
string sequencer(int size, int seedNum);

int main(){
    Tester tester;
    VDetect* vdetect;
    bool result;

    cout << "**********************" << endl;
    cout << "***** BEGIN TEST *****" << endl;
    cout << "**********************" << endl << endl;

    {
        cout << "\n*** TEST BLOCK ONE***" << endl << endl;
        cout << "This will test normal, noncolliding insertion to ensure it works properly" << endl << endl;
        vdetect = new VDetect(MINPRIME, hashCode, DOUBLEHASH);

        cout << "Testing Non-colliding Insertions:\n\t";
        result = tester.testInsertion(*vdetect);
        tester.testCondition(result);


        delete vdetect;
        cout << "\n***END TEST BLOCK ONE***" << endl;  
    }
    {
        cout << "\n*** TEST BLOCK TWO***" << endl << endl;
        cout << "This will test collision cases of insertion, to see if it is handled properly" << endl << endl;
        
        //first test is on NONE probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, NONE);

        cout << "Testing Colloisions with the NONE Policy:\n\t";
        result = tester.collisionInsertion(*vdetect);
        tester.testCondition(result);
        delete vdetect;

        //second test is on DOUBLEHASH probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, DOUBLEHASH);

        cout << "Testing Colloisions with the DOUBLEHASH Policy:\n\t";
        result = tester.collisionInsertion(*vdetect);
        tester.testCondition(result);
        delete vdetect; 

        //final test is with the quadratic policy
        //second test is on DOUBLEHASH probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, QUADRATIC);

        cout << "Testing Colloisions with the QUADRATIC Policy:\n\t";
        result = tester.collisionInsertion(*vdetect);
        tester.testCondition(result);
        delete vdetect;         

        cout << "\n***END TEST BLOCK TWO ***" << endl;         
    }
    {
        cout << "\n*** TEST BLOCK THREE***" << endl << endl;
        cout << "This will test collision cases of getVirus, to see if it is handled properly" << endl << endl;
        
        //first test is on NONE probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, NONE);

        cout << "Testing getVirus Colloisions with the NONE Policy:\n\t";
        result = tester.collisionGetVirus(*vdetect);
        tester.testCondition(result);
        delete vdetect;
        
        //second test is on DOUBLEHASH probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, DOUBLEHASH);

        cout << "Testing getVirus Collisions with the DOUBLEHASH Policy:\n\t";
        result = tester.collisionGetVirus(*vdetect);
        tester.testCondition(result);
        delete vdetect; 
    
        //final test is with the quadratic policy
        vdetect = new VDetect(MINPRIME, hashCode, QUADRATIC);

        cout << "Testing getVirus Collisions with the QUADRATIC Policy:\n\t";
        result = tester.collisionGetVirus(*vdetect);
        tester.testCondition(result);
        delete vdetect;         

        cout << "\n***END TEST BLOCK THREE***" << endl;
    }
    {
        cout << "\n*** TEST BLOCK FOUR***" << endl << endl;
        cout << "This will test normal, noncolliding getVirus to ensure it works properly" << endl << endl;
        vdetect = new VDetect(MINPRIME, hashCode, DOUBLEHASH);

        cout << "Testing Non-colliding getVirus:\n\t";
        result = tester.normalGetVirus(*vdetect);
        tester.testCondition(result);

        delete vdetect;
        cout << "\n***END TEST BLOCK FOUR***" << endl;         
    }
    {
        cout << "\n*** TEST BLOCK FIVE***" << endl << endl;
        cout << "This will test error cases of getVirus to ensure it recognizes error" << endl << endl;
        vdetect = new VDetect(MINPRIME, hashCode, NONE);

        cout << "Testing error cases of getVirus:\n\t";
        result = tester.errorGetVirus(*vdetect);
        tester.testCondition(result);

        delete vdetect;
        cout << "\n***END TEST BLOCK FOUR***" << endl;            
    }
    {
        cout << "\n*** TEST BLOCK SIX***" << endl << endl;
        cout << "This will test removal of non-colliding keys" << endl << endl;
        vdetect = new VDetect(MINPRIME, hashCode, NONE);

        cout << "Testing non collisions of remove:\n\t";
        result = tester.normalRemove(*vdetect);
        tester.testCondition(result);

        delete vdetect;
        cout << "\n***END TEST BLOCK SIX***" << endl;            
    }
    {
        cout << "\n*** TEST BLOCK SEVEN***" << endl << endl;
        cout << "This will test collision cases of removal, to see if it is handled properly" << endl << endl;
        
        //first test is on NONE probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, NONE);

        cout << "Testing removal with the NONE Policy:\n\t";
        result = tester.collisionRemove(*vdetect);
        tester.testCondition(result);
        delete vdetect;
        
        //second test is on DOUBLEHASH probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, DOUBLEHASH);

        cout << "Testing removal with the DOUBLEHASH Policy:\n\t";
        result = tester.collisionRemove(*vdetect);
        tester.testCondition(result);
        delete vdetect; 

        //final test is with the quadratic policy
        //second test is on DOUBLEHASH probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, QUADRATIC);

        cout << "Testing removal with the QUADRATIC Policy:\n\t";
        result = tester.collisionRemove(*vdetect);
        tester.testCondition(result);
        delete vdetect;         

        cout << "\n***END TEST BLOCK SEVEN ***" << endl;
    }
    {
        cout << "\n*** TEST BLOCK EIGHT***" << endl << endl;
        cout << "This will test to see if rehashing is triggered properly after insertion" << endl << endl;
        
        //first test is on NONE probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, NONE);

        cout << "Testing insertRehash with the NONE Policy:\n\t";
        result = tester.insertRehash(*vdetect);
        tester.testCondition(result);
        delete vdetect;
        
        //second test is on DOUBLEHASH probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, DOUBLEHASH);

        cout << "Testing insertRehash with the DOUBLEHASH Policy:\n\t";
        result = tester.insertRehash(*vdetect);
        tester.testCondition(result);
        delete vdetect; 

        //final test is with the quadratic policy
        //second test is on DOUBLEHASH probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, QUADRATIC);

        cout << "Testing insertRehash with the QUADRATIC Policy:\n\t";
        result = tester.insertRehash(*vdetect);
        tester.testCondition(result);
        delete vdetect;         

        cout << "\n***END TEST BLOCK EIGHT ***" << endl;        
    }
    {
        cout << "\n*** TEST BLOCK NINE***" << endl << endl;
        cout << "This will test to see if that everything is properly reordered after being rehashed" << endl << endl;
        
        //first test is on NONE probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, NONE);

        cout << "Testing completeLambdaRehash with the NONE Policy:\n\t";
        result = tester.completeLambdaRehash(*vdetect);
        tester.testCondition(result);
        delete vdetect;
        
        //second test is on DOUBLEHASH probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, DOUBLEHASH);

        cout << "Testing completeLambdaRehash with the DOUBLEHASH Policy:\n\t";
        result = tester.completeLambdaRehash(*vdetect);
        tester.testCondition(result);
        delete vdetect; 

        //final test is with the quadratic policy
        //second test is on DOUBLEHASH probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, QUADRATIC);

        cout << "Testing insertRehash with the QUADRATIC Policy:\n\t";
        result = tester.completeLambdaRehash(*vdetect);
        tester.testCondition(result);
        delete vdetect;         

        cout << "\n***END TEST BLOCK NINE ***" << endl;               
    }
    {
        cout << "\n*** TEST BLOCK TEN***" << endl << endl;
        cout << "This will test to see if rehashing triggers after numerous removals" << endl << endl;
        
        //first test is on NONE probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, NONE);

        cout << "Testing removeRehash with the NONE Policy:\n\t";
        result = tester.removeRehash(*vdetect);
        tester.testCondition(result);
        delete vdetect;
        
        //second test is on DOUBLEHASH probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, DOUBLEHASH);

        cout << "Testing removeRehash with the DOUBLEHASH Policy:\n\t";
        result = tester.removeRehash(*vdetect);
        tester.testCondition(result);
        delete vdetect; 

        //final test is with the quadratic policy
        //second test is on DOUBLEHASH probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, QUADRATIC);

        cout << "Testing removeRehash with the QUADRATIC Policy:\n\t";
        result = tester.removeRehash(*vdetect);
        tester.testCondition(result);
        delete vdetect;         

        cout << "\n***END TEST BLOCK TEN ***" << endl;         
    }
    {
        cout << "\n*** TEST BLOCK ELEVEN***" << endl << endl;
        cout << "This will test to see if rehashing triggers after numerous removals" << endl << endl;
        
        //first test is on NONE probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, NONE);

        cout << "Testing completeDeleteRehash with the NONE Policy:\n\t";
        result = tester.completeDeleteRehash(*vdetect);
        tester.testCondition(result);
        delete vdetect;
        
        //second test is on DOUBLEHASH probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, DOUBLEHASH);

        cout << "Testing completeDeleteRehash with the DOUBLEHASH Policy:\n\t";
        result = tester.completeDeleteRehash(*vdetect);
        tester.testCondition(result);
        delete vdetect; 

        //final test is with the quadratic policy
        //second test is on DOUBLEHASH probing, so repeats are rid of
        vdetect = new VDetect(MINPRIME, hashCode, QUADRATIC);

        cout << "Testing completeDeleteRehash with the QUADRATIC Policy:\n\t";
        result = tester.completeDeleteRehash(*vdetect);
        tester.testCondition(result);
        delete vdetect;         

        cout << "\n***END TEST BLOCK ELEVEN***" << endl;
    }


    cout << "\n**********************" << endl;
    cout << "***** END TEST *******" << endl;
    cout << "**********************" << endl;

    return 0;
}

unsigned int hashCode(const string str) {
   unsigned int val = 0 ;
   const unsigned int thirtyThree = 33 ;  // magic number from textbook
   for (int i = 0 ; i < int(str.length()); i++)
      val = val * thirtyThree + str[i] ;
   return val ;
}

string sequencer(int size, int seedNum){
    //this function returns a random DNA sequence
    string sequence = "";
    Random rndObject(0,3);
    rndObject.setSeed(seedNum);
    for (int i=0;i<size;i++){
        sequence = sequence + ALPHA[rndObject.getRandNum()];
    }
    return sequence;
}

//testCondition
//displays the output of the test whether it passed or failed
void Tester::testCondition(bool var){
    if (var){
        cout << "This test has passed" << endl << endl;
    }
    else{
        cout << "This test has failed" << endl << endl;
    }
}

//testInsertion
//tests normal insertion with no collisions
bool Tester::testInsertion(VDetect& vdetect){
    bool result = true;
    Random RndID(MINID,MAXID);
    int trueSize = 0;

    //for loop goes through and adds 49 virus objects, which is enough to not trigger the load factor
    for (int i = 0; i < 49; i++){

        //virus created with a random id, the sequencer creates a random DNA sequence
        Virus virus = Virus(sequencer(5,i), RndID.getRandNum());

        //first the number hashed to is checked, and the potential index number is checked

        int potentialIndex = hashCode(virus.getKey()) % vdetect.m_currentCap;

        //if this index is taken already, then it is skipped as not to mess up the results of this test
        if ((vdetect.m_currentTable[potentialIndex] == EMPTY)){
            //else, if it is not, then variables are checked
            int oldSize = vdetect.m_currentSize;

            //the virus is inserted
            result = result && vdetect.insert(virus);

            //it checks to ensure member variables are changed
            result = result && (vdetect.m_currentSize == (oldSize + 1));
            result = result && (vdetect.m_currentTable[potentialIndex] == virus);

            //it is then checked to ensure it is found
            result = result && (vdetect.getVirus(virus.getKey(), virus.getID()) == virus);

            ++trueSize; //trueSize increased by one
        }
    }

    //result checks to ensure the size is correct
    result = result && (vdetect.m_currentSize == trueSize);

    return result; //return result
}

//collisionInsertion
//tests inserting colliding keys
bool Tester::collisionInsertion(VDetect& vdetect){
    //variables which are used throughout the testing process
    bool result = true;
    int oldSize;
    int potentialIndex;
    int trueSize = 1;

    //newVirus created, will be used for hashing
    Virus* newVirus = new Virus("AGG", MINID);
    
    //hashNumber is saved
    int hashNum = hashCode(newVirus->getKey()) % vdetect.m_currentCap;
    
    //the object is then inserted
    result = result && vdetect.insert(*newVirus);

    //for loop goes through and checks to ensure that collisions are handled properly depending on which case is used
    for (int i = 1; i < 49; i++){
        oldSize = vdetect.m_currentSize;
        newVirus->setID(MINID + i); //ID updated so that ID is always unique

        //first checked if it is none,
        //if it is none, then nothing should change
        if (vdetect.m_currProbing == NONE){
            result = result && !(vdetect.insert(*newVirus));
            result = result && (vdetect.getVirus(newVirus->getKey(), newVirus->getID()) == EMPTY);
            result = result && (vdetect.m_currentSize == oldSize);
        }
        //else, it checks if double probing works properly
        else if (vdetect.m_currProbing == DOUBLEHASH){
            potentialIndex = ((hashNum + (i * (11-(vdetect.m_hash(newVirus->getKey()) % 11)))) % vdetect.m_currentCap);
            
            //if statement checks to ensure that the potential index spot is empty
            if (vdetect.m_currentTable[potentialIndex] == EMPTY){
                result = result & vdetect.insert(*newVirus);

                //it checks to ensure member variables are changed
                result = result && (vdetect.m_currentSize == (oldSize + 1));
                result = result && (vdetect.m_currentTable[potentialIndex] == *newVirus);

                //it is then checked to ensure it is found
                result = result && (vdetect.getVirus(newVirus->getKey(), newVirus->getID()) == *newVirus);
        
                ++trueSize; //trueSize increased by one
            }
        }
        //else, it is quadratic probing, so it's checked to ensure it is right
        else{
            potentialIndex = (hashNum + (i * i)) % vdetect.m_currentCap;
            //if statement checks to ensure that the potential index spot is empty
            if (vdetect.m_currentTable[potentialIndex] == EMPTY){
                result = result & vdetect.insert(*newVirus);
                //it checks to ensure member variables are changed
                result = result && (vdetect.m_currentSize == (oldSize + 1));
                result = result && (vdetect.m_currentTable[potentialIndex] == *newVirus);

                //it is then checked to ensure it is found
                result = result && (vdetect.getVirus(newVirus->getKey(), newVirus->getID()) == *newVirus);
                ++trueSize; //trueSize increased by one
            }
        }
    }

    result = result && (vdetect.m_currentSize == trueSize); //result checks to ensure that currentSize is updated

    delete newVirus; //newVirus deleted
    return result;
}

//collisionGetVirus
//checks to ensure collisions are handled properly with getVirus
bool Tester::collisionGetVirus(VDetect& vdetect){
    bool result = true;
    Virus* virus = new Virus("AGG", MINID);
    //first thing done is to add insertCollision collisions
    result = result && collisionInsertion(vdetect);

    //next, getVirus is tested for each of them
    //THIS TESTS BEFORE A REHASH IS TRIGGERED
    for (int i = 0; i < 49; i++){
        virus->setID(MINID + i);
        result = result && (vdetect.getVirus(virus->getKey(),virus->getID()) == *virus);
        //if currProbing is NONE, breaks early
        if (vdetect.m_currProbing == NONE){
            break;
        }
    }


    //triggers a rehash specifically for quadratic, since otherwise the limit will be hit
    if (vdetect.m_currProbing == QUADRATIC){
        virus->setKey("TTA");
        result = result && vdetect.insert(*virus);
        //checks to ensure that size is correct and that it was inserted
        //getVirus tested each time
        result = result && (vdetect.getVirus(virus->getKey(),virus->getID()) == *virus);
        virus->setKey("AGG");
    }

    //THIS TESTS GETVIRUS AFTER A REHASH IS EXPLICITLY TRIGGERED
    if (vdetect.m_currProbing != NONE){
        //next, more insertions are done, such that they end up in the old table
        //101 is used so that a rehash isn't triggered
        for (int i = 49; i < 101; i++){
            virus->setID(MINID + i);
            result = result && vdetect.insert(*virus);
            //checks to ensure that size is correct and that it was inserted
            //getVirus tested each time
            result = result && (vdetect.getVirus(virus->getKey(),virus->getID()) == *virus);
        }
        
        //next, getVirus is tested for each of them to ensure everything exists
        for (int i = 0; i < 101; i++){
            virus->setID(MINID + i);
            result = result && (vdetect.getVirus(virus->getKey(),virus->getID()) == *virus);
            //if currProbing is NONE, breaks early
            if (vdetect.m_currProbing == NONE){
                break;
            }
        }
    }

    delete virus;
    return result;
}

//normalGetVirus
//tests the hash table's getVirus for non-colliding keys
bool Tester::normalGetVirus(VDetect& vdetect){
    bool result = true;
    string compare[49]; //array string used to compare keys
    //noncolliding keys are inserted, and then are tested by getVirus to ensure they exist
    int trueSize = 0;

    Virus* virus = new Virus("G", MINID);

    //for loop goes through and adds 49 virus objects to see if getVirus will work for everything
    //49 is used to prevent rehashing, but repeats will occur
    for (int i = 0; i < 49; i++){
        //key and setID changed each time to ensure no collisions occur
        virus->setKey(sequencer(5, i));
        virus->setID(MINID + i);
        
        compare[i] = virus->getKey();
         
        //first the number hashed to is checked, and the potential index number is checked
        int potentialIndex = hashCode(virus->getKey()) % vdetect.m_currentCap;

        //if this index is taken already, then it is skipped as not to mess up the results of this test
        if ((vdetect.m_currentTable[potentialIndex] == EMPTY)){

            //else, if it is not, then variables are checked
            int oldSize = vdetect.m_currentSize;
            //the virus is inserted
            result = result && vdetect.insert(*virus);

            //it checks to ensure member variables are changed
            result = result && ((vdetect.m_currentSize + vdetect.m_oldSize) == (oldSize + 1));
            result = result && (vdetect.m_currentTable[potentialIndex] == *virus);

            //it is then checked to ensure it is found
            result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == *virus);

            ++trueSize; //trueSize increased by one
        }
        else{
            compare[i] = "EMPTY"; //else, compare i set to empty, to make it easier to use getVirus
        }
    }
    //result checks to ensure the size is correct
    result = result && ((vdetect.m_currentSize + vdetect.m_oldSize) == trueSize);

    //finally, getVirus is checked again to ensure it works properly
    for (int i = 0; i < 20; i++){
        virus->setKey(compare[i]);
        virus->setID(MINID + i);
        if (virus->getKey() != "EMPTY"){
            result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == *virus);
        }
    }

    delete virus;
    return result;
}

//errorGetVirus
//tests getVirus cases for errors
bool Tester::errorGetVirus(VDetect& vdetect){
    bool result = true;
    Virus* virus = new Virus("AGG", MINID);

    //virus is first inserted into the hash table
    result = result && vdetect.insert(*virus);

    //it is inserted again, and since it is a NONE hashing policy, it should not be inserted nor found the second time
    virus->setID(MAXID);
    result = result && !(vdetect.insert(*virus));
    result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == EMPTY);

    //next test done is attempting to insert a virus with an ineligible ID
    virus->setKey("ATT");
    virus->setID(-1);
    result = result && !(vdetect.insert(*virus));
    result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == EMPTY);

    //it is then tested again with an ID which is too large
    virus->setID(10000000);
    result = result && !(vdetect.insert(*virus));
    result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == EMPTY);

    //next, we try inserting a key which is too long to be inserted
    virus->setID(MAXID);
    virus->setKey("AGGATA");
    result = result && !(vdetect.insert(*virus));
    result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == EMPTY);

    //next getVirus test is a valid Virus which was never inserted
    virus->setID(MAXID);
    virus->setKey("CTG");
    result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == EMPTY);
    
    //finally, a virus is added, then deleted
    result = result && vdetect.insert(*virus);
    result = result && vdetect.remove(*virus);
    result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == EMPTY);

    delete virus;
    return result;
}

//normal remove
//tests removal for non colliding keys
bool Tester::normalRemove(VDetect& vdetect){
    bool result = true;
    string compare[49]; //array string used to compare keys
    //noncolliding keys are inserted, and then are tested by getVirus to ensure they exist
    Virus* virus = new Virus("G", MINID);
    int oldDeleted;
    int oldSize;

    //for loop goes through 49 times
    //20 times because the seeding triggers a repeat at i = 49, so easier to avoid it
    for (int i = 0; i < 49; i++){
        //key and setID changed each time to ensure no collisions occur
        virus->setKey(sequencer(5, i));
        virus->setID(MINID + i);
        
        compare[i] = virus->getKey();
         
        //first the number hashed to is checked, and the potential index number is checked
        int potentialIndex = hashCode(virus->getKey()) % vdetect.m_currentCap;

        //if this index is taken already, then it is skipped as not to mess up the results of this test
        if ((vdetect.m_currentTable[potentialIndex] == EMPTY)){

            //else, if it is not, then variables are checked
            oldSize = vdetect.m_currentSize;
            //the virus is inserted
            result = result && vdetect.insert(*virus);

            //it checks to ensure member variables are changed
            result = result && ((vdetect.m_currentSize + vdetect.m_oldSize) == (oldSize + 1));
            result = result && (vdetect.m_currentTable[potentialIndex] == *virus);

            //it is then checked to ensure it is found
            result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == *virus);
        }
        else{
            compare[i] = "EMPTY"; //else, compare i set to empty, to make it easier to remove
        }
    }

    int bound = 0;
    int numDeleted = 0;

    for (int i = 0; i < 49; i++){
        if (compare[i] != "EMPTY"){
            ++bound;
        }
    }
    //next, removal is tested, to ensure that it works properly
    //only thirty are attempted to be removed, as to not trigger the deleted ratio
    for (int i = 0; i < 30; i++){
        virus->setKey(compare[i]);
        virus->setID(MINID + i);
        oldSize = vdetect.m_currentSize;
        oldDeleted = vdetect.m_currNumDeleted;
        
        if (compare[i] != "EMPTY"){
            result = result && vdetect.remove(*virus);

            //checks to ensure member variables updated
            result = result && (vdetect.m_currNumDeleted == (oldDeleted + 1));
            //checks to ensure getVirus cannot find the virus in question
            result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == EMPTY);
            ++numDeleted;
        }
    }
    //checks to ensure all the member variables are properly set
    result = result && (vdetect.m_currentSize == bound);
    result = result && (vdetect.m_currNumDeleted == numDeleted);

    //checks to ensure the deleted ratio is properly set
    result = result && (vdetect.deletedRatio() == float(numDeleted)/float(bound));
    
    delete virus;
    return result;
}

//collisionRemove
//tests remove to see if its accurate
bool Tester::collisionRemove(VDetect& vdetect){
    bool result = true;
    Virus* virus = new Virus("AGG", MINID);
    int oldDeleted;
    int numDeleted = 0;

    result = result && collisionInsertion(vdetect); //collisionInsertion vdetect called to fill the tables
    
    //if currProbing is none, another one is inserted to avoid rehashing
    if (vdetect.m_currProbing == NONE){
        virus->setKey("CAAG");
        result = result && vdetect.insert(*virus);
        virus->setKey("AGG"); //returned back to AGG
    }

    //for loop goes through and tries to remove everything currently in each table
    //up to 39 to prevent a rehash

    for (int i = 0; i < 39; i++){
        virus->setID(MINID + i);
        oldDeleted = vdetect.m_currNumDeleted;
        
        result = result && vdetect.remove(*virus); //virus is attempted to be removed

       //checks to ensure member variables updated
        result = result && ((vdetect.m_currNumDeleted + vdetect.m_oldNumDeleted) == (oldDeleted + 1));

        //checks to ensure getVirus cannot find the virus in question
        result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == EMPTY);
        ++numDeleted;

        if (vdetect.m_currProbing == NONE){
            break;
        }
    }

    //if statement checks which currProbing policy is used
    //if NONE, checks that size is 1 since only 1 removed, while opposite is true for QUADRATIC/DOUBLE
    if (vdetect.m_currProbing == NONE){
        //checks to ensure all member variables are properly set
        result = result && (vdetect.m_currentSize == 2);
        result = result && (vdetect.m_currNumDeleted == 1);
        result = result && (vdetect.deletedRatio() == 0.5); //should be 0.5, since 1/2 = 0
    }
    else{
        result = result && (vdetect.m_currentSize == 49);
        result = result && (vdetect.m_currNumDeleted == 39);
        result = result && (vdetect.deletedRatio() == float(39)/float(49));
    }

    delete virus;
    return result;
}

//insertRehash
//Attempts to see if a rehash is triggered once a certain number of insertions are reached
bool Tester::insertRehash(VDetect& vdetect){
    bool result = true;
    Virus* virus = new Virus("", MINID);
    Virus checkForAll[51];


    for (int i = 0; i < 1000; i++){
        //virus created with a random id, the sequencer creates a random DNA sequence
        virus->setKey(sequencer(5, i));
        virus->setID(MINID + i);

        //first the number hashed to is checked, and the potential index number is checked
        int potentialIndex = hashCode(virus->getKey()) % vdetect.m_currentCap;

        //if this index is taken already, then it is skipped as not to mess up the results of this test
        if ((vdetect.m_currentTable[potentialIndex] == EMPTY)){
            //the virus is inserted
            result = result && vdetect.insert(*virus);
            
            //it is then checked to ensure it is found
            result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == *virus);
            checkForAll[vdetect.m_currentSize - 1] = *virus;
        }
        
        if (vdetect.m_currentSize == 50){
            break;
        }
    }

    vdetect.changeProbPolicy(QUADRATIC); //always changes to a qaudratic prob policy

    //next, we are going to try triggering the rehash to ensure it triggers properly
    //high number used specifically for NONE
    for (int i = 100; i < MAXID; i++){
        //virus created with a random id, the sequencer creates a random DNA sequence
        virus->setKey(sequencer(5, i));
        virus->setID(MINID + i);

        //first the number hashed to is checked, and the potential index number is checked
        int potentialIndex = hashCode(virus->getKey()) % vdetect.m_currentCap;

        //if this index is taken already, then it is skipped as not to mess up the results of this test
        if ((vdetect.m_currentTable[potentialIndex] == EMPTY)){
            int newCheck = vdetect.findNextPrime(4 *((vdetect.m_currentSize+1) - vdetect.m_currNumDeleted));
            //the virus is inserted
            result = result && vdetect.insert(*virus);

            //next, we check that oldTable is no longer a nullptr, and that everything transferred over
            result = result && (vdetect.m_oldTable != nullptr);
            result = result && (vdetect.m_oldSize == 51);
            result = result && (vdetect.m_currentSize == 12);
            result = result && (vdetect.m_oldNumDeleted == 12);
            result = result && (vdetect.m_currProbing == QUADRATIC);
            result = result && (vdetect.lambda() < 0.5); //checks to ensure lambda was reset
            result = result && (vdetect.m_currentCap == newCheck);
            //it is then checked to ensure it is found
            result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == *virus);
            checkForAll[50] = *virus;
            break;
        }
    }


    //checks to ensure everything exists within the list
    for (int i = 0; i < 51; i++){
        result = result && (vdetect.getVirus(checkForAll[i].getKey(), checkForAll[i].getID()) == checkForAll[i]);
    }
    delete virus;
    return result;
}

//completeLambdaRehash
//Tests the lamba after a large number of insertion
bool Tester::completeLambdaRehash(VDetect& vdetect){
    bool result = true;
    Virus* virus = new Virus("", MINID);
    Virus checkForAll[55];

    //first for loop goes through and adds the maximum before it is triggered
    for (int i = 0; i < 1000; i++){
        //virus created with a random id, the sequencer creates a random DNA sequence
        virus->setKey(sequencer(5, i));
        virus->setID(MINID + i);

        //first the number hashed to is checked, and the potential index number is checked
        int potentialIndex = hashCode(virus->getKey()) % vdetect.m_currentCap;

        //if this index is taken already, then it is skipped as not to mess up the results of this test
        if ((vdetect.m_currentTable[potentialIndex] == EMPTY)){
            //the virus is inserted
            result = result && vdetect.insert(*virus);
            
            //it is then checked to ensure it is found
            result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == *virus);
            checkForAll[vdetect.m_currentSize - 1] = *virus;
        }
        
        if (vdetect.m_currentSize == 50){
            break;
        }
    }

    vdetect.changeProbPolicy(QUADRATIC); //always changes to a qaudratic prob policy

    int count = 0; //count keeps track of how many times we've gone through the for loop

    //next, we are going to try triggering the rehash to ensure it triggers properly
    //high number used specifically for NONE
    for (int i = 100; i < MAXID; i++){
        //virus created with a random id, the sequencer creates a random DNA sequence
        virus->setKey(sequencer(5, i));
        virus->setID(MINID + i);

        //first the number hashed to is checked, and the potential index number is checked
        int potentialIndex = hashCode(virus->getKey()) % vdetect.m_currentCap;

        //if this index is taken already, then it is skipped as not to mess up the results of this test
        if ((vdetect.m_currentTable[potentialIndex] == EMPTY)){
            //the virus is inserted
            result = result && vdetect.insert(*virus);
    
            //it is then checked to ensure it is found
            result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == *virus);
            checkForAll[(50 + count)] = *virus;

            ++count;
            //once 5 take place, everything should be reset, so that is checked
            if (count == 5){
                result = result && (vdetect.m_oldTable == nullptr);
                result = result && (vdetect.m_oldSize == 0);
                result = result && (vdetect.m_oldNumDeleted == 0); 
                result = result && (vdetect.m_oldCap == 0);
                result = result && (vdetect.m_oldProbing == NONE);
 
                //checks to ensure lambda is still less than .5 since everything has been transferred over by this point
                result = result && (vdetect.lambda() < 0.5); //checks to ensure lambda was reset
                break; //breaks
            }
        }
    }

    //checks to ensure everything exists within the list
    for (int i = 0; i < 55; i++){
        result = result && (vdetect.getVirus(checkForAll[i].getKey(), checkForAll[i].getID()) == checkForAll[i]);
    }
    delete virus;
    return result;
}

//removeRehash
//tests to see if rehash works properly after the removal threshold is reached
bool Tester::removeRehash(VDetect& vdetect){
    bool result = true;
    Virus* virus = new Virus("", MINID);
    Virus checkForAll[51];

    //for loop goes through and add viruses
    for (int i = 0; i < 1000; i++){
        //virus created with a random id, the sequencer creates a random DNA sequence
        virus->setKey(sequencer(5, i));
        virus->setID(MINID + i);

        //first the number hashed to is checked, and the potential index number is checked
        int potentialIndex = hashCode(virus->getKey()) % vdetect.m_currentCap;

        //if this index is taken already, then it is skipped as not to mess up the results of this test
        if ((vdetect.m_currentTable[potentialIndex] == EMPTY)){
            //the virus is inserted
            result = result && vdetect.insert(*virus);
            
            //it is then checked to ensure it is found
            result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == *virus);
            checkForAll[vdetect.m_currentSize - 1] = *virus;
        }
        
        if (vdetect.m_currentSize == 50){
            break;
        }
    }

    vdetect.changeProbPolicy(QUADRATIC); //always changes to a qaudratic prob policy

    //for loop goes through and tries to remove everything currently in each table
    //breaks when the deleted ratio is triggered
    for (int i = 0; i < 50; i++){
        result = result && vdetect.remove(checkForAll[i]); //virus is attempted to be removed
        
        //checks to ensure getVirus cannot find the virus in question
        result = result && (vdetect.getVirus(checkForAll[i].getKey(), checkForAll[i].getID()) == EMPTY);

        //if i equals 37, then we insert again to trigger a rehash
        if (i == 38){            
            int newCheck = vdetect.findNextPrime(4 *((vdetect.m_currentSize) - (vdetect.m_currNumDeleted + 1)));

            result = result && vdetect.remove(checkForAll[i + 1]); //virus is attempted to be removed

            //next, we check that oldTable is no longer a nullptr, and that everything transferred over
            result = result && (vdetect.m_oldTable != nullptr);

            //if statements check to ensure member variables are the correct values
            result = result && (vdetect.m_oldSize == 50);
            result = result && (vdetect.m_oldNumDeleted == 40 + vdetect.m_currentSize);
            result = result && (vdetect.m_currentSize == 2);
            result = result && (vdetect.m_currNumDeleted == 0);
            result = result && (vdetect.m_currProbing == QUADRATIC);
            result = result && (vdetect.m_currentCap == newCheck);

            //checks to see if the deleted ratio factor is now zero
            result = result && (vdetect.deletedRatio() == 0); 
            break;
        }
    }

    //checks to ensure everything exists within the list
    for (int i = 0; i < 50; i++){
        if (i < 40){
            result = result && (vdetect.getVirus(checkForAll[i].getKey(), checkForAll[i].getID()) == EMPTY);
        }
        else{
            result = result && (vdetect.getVirus(checkForAll[i].getKey(), checkForAll[i].getID()) == checkForAll[i]);
        } 
    }

    delete virus;
    return result;
}

//completeDeleteRehash
//Tests complete deletion rehash after triggered by deleted ratio
bool Tester::completeDeleteRehash(VDetect& vdetect){
    bool result = true;
    Virus* virus = new Virus("", MINID);
    Virus checkForAll[50];

    //for loop goes through and add viruses
    for (int i = 0; i < 1000; i++){
        //virus created with a random id, the sequencer creates a random DNA sequence
        virus->setKey(sequencer(5, i));
        virus->setID(MINID + i);

        //first the number hashed to is checked, and the potential index number is checked
        int potentialIndex = hashCode(virus->getKey()) % vdetect.m_currentCap;

        //if this index is taken already, then it is skipped as not to mess up the results of this test
        if ((vdetect.m_currentTable[potentialIndex] == EMPTY)){
            //the virus is inserted
            result = result && vdetect.insert(*virus);
            
            //it is then checked to ensure it is found
            result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == *virus);
            checkForAll[vdetect.m_currentSize - 1] = *virus;
        }
        
        if (vdetect.m_currentSize == 50){
            break;
        }
    }

    vdetect.changeProbPolicy(QUADRATIC); //always changes to a qaudratic prob policy

    //for loop goes through and tries to remove everything currently in each table
    //breaks when the deleted ratio is triggered
    for (int i = 0; i < 40; i++){
        result = result && vdetect.remove(checkForAll[i]); //virus is attempted to be removed
        
        //checks to ensure getVirus cannot find the virus in question
        result = result && (vdetect.getVirus(checkForAll[i].getKey(), checkForAll[i].getID()) == EMPTY);
        checkForAll[i] = EMPTY; //checkForAll set to empty
        
        //if i equals 39, we check to ensure that a rehash was triggered
        if (i == 39){
            //next, we check that oldTable is no longer a nullptr, and that everything transferred over
            result = result && (vdetect.m_oldTable != nullptr);

            //if statements check to ensure member variables are the correct values
            result = result && (vdetect.m_oldSize == 50);
            result = result && (vdetect.m_oldNumDeleted == 40 + vdetect.m_currentSize);
            result = result && (vdetect.m_currentSize == 2);
            result = result && (vdetect.m_currNumDeleted == 0);
            result = result && (vdetect.m_currProbing == QUADRATIC);

            //checks to see if the deleted ratio factor is now zero
            result = result && (vdetect.deletedRatio() == 0);
        }
    }

    int count = 0;

    //first we try inserting four times to ensure that the rehash empties the table
    for (int i = 0; i < 1000; i++){
        //virus created with a random id, the sequencer creates a random DNA sequence
        virus->setKey(sequencer(5, i));
        virus->setID(MINID + i);

        //first the number hashed to is checked, and the potential index number is checked
        int potentialIndex = hashCode(virus->getKey()) % vdetect.m_currentCap;

        //if this index is taken already, then it is skipped as not to mess up the results of this test
        if ((vdetect.m_currentTable[potentialIndex] == EMPTY)){
            result = result && (vdetect.m_oldTable != nullptr);
            //the virus is inserted
            result = result && vdetect.insert(*virus);
            
            //it is then checked to ensure it is found
            result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == *virus);
            checkForAll[count] = *virus;
            ++count;

            //if count equals 4, everything should then be cleared out by this point
            if (count == 4){
                result = result && (vdetect.m_oldTable == nullptr);
                result = result && (vdetect.m_oldSize == 0);
                result = result && (vdetect.m_oldNumDeleted == 0); 
                result = result && (vdetect.m_oldCap == 0);
                result = result && (vdetect.m_oldProbing == NONE);
                break; //breaks
            }
        }
    }

    count = 0; //count reset to zero

    //next, we will insert, then remove again, then trigger the deleted rehash through removal
    //for loop goes through and add viruses
    for (int i = 0; i < 1000; i++){
        //virus created with a random id, the sequencer creates a random DNA sequence
        virus->setKey(sequencer(5, i));
        virus->setID(MINID + i);

        //first the number hashed to is checked, and the potential index number is checked
        int potentialIndex = hashCode(virus->getKey()) % vdetect.m_currentCap;

        //if this index is taken already, then it is skipped as not to mess up the results of this test
        if ((vdetect.m_currentTable[potentialIndex] == EMPTY)){
            //the virus is inserted
            result = result && vdetect.insert(*virus);
            ++count;
            //it is then checked to ensure it is found
            result = result && (vdetect.getVirus(virus->getKey(), virus->getID()) == *virus);
            checkForAll[3 + count] = *virus;
        }
        
        if (vdetect.m_currentSize == 50){
            break;
        }
    }

    //next, we attempt to remove and trigger another rehash
    vdetect.changeProbPolicy(DOUBLEHASH); //always changes to a qaudratic prob policy
    result = result && (vdetect.m_currentSize == 50);

    //now, we will remove again, and trigger a rehash
    for (int i = 0; i < 40; i++){
        result = result && vdetect.remove(checkForAll[i]); //virus is attempted to be removed
        
        //checks to ensure getVirus cannot find the virus in question
        result = result && (vdetect.getVirus(checkForAll[i].getKey(), checkForAll[i].getID()) == EMPTY);
        checkForAll[i] = EMPTY; //checkForAll set to empty
        
        //if i equals 39, we check to ensure that a rehash was triggered
        if (i == 39){
            //next, we check that oldTable is no longer a nullptr, and that everything transferred over
            result = result && (vdetect.m_oldTable != nullptr);

            //if statements check to ensure member variables are the correct values
            result = result && (vdetect.m_oldSize == 50);
            result = result && (vdetect.m_oldNumDeleted == 40 + vdetect.m_currentSize);
            result = result && (vdetect.m_currentSize == 2);
            result = result && (vdetect.m_currNumDeleted == 0);
            result = result && (vdetect.m_currProbing == DOUBLEHASH);

            //checks to see if the deleted ratio factor is now zero
            result = result && (vdetect.deletedRatio() == 0);
        }
    }

    //next, we will remove again to check and ensure that the rehash is triggered
    for (int i = 40; i < 43; i++){
        result = result && (vdetect.m_oldTable != nullptr); //checks to ensure our table is not nullptr
        result = result && vdetect.remove(checkForAll[i]); //virus is attempted to be removed
        //checks to ensure getVirus cannot find the virus in question
        result = result && (vdetect.getVirus(checkForAll[i].getKey(), checkForAll[i].getID()) == EMPTY);
        checkForAll[i] = EMPTY;
    }  

    //we then check to ensure the member variables are accurate
    result = result && (vdetect.m_oldTable == nullptr);
    result = result && (vdetect.m_oldSize == 0);
    result = result && (vdetect.m_oldNumDeleted == 0); 
    result = result && (vdetect.m_oldCap == 0);
    result = result && (vdetect.m_oldProbing == NONE);

    //finally, we check to see that only them proper checkForAll array members are there
    for (int i = 0; i < 50; i++){
        //cout << i << " " << result << endl;
        if (i < 43){
            result = result && (vdetect.getVirus(checkForAll[i].getKey(), checkForAll[i].getID()) == EMPTY);
        }
        else{
            result = result && (vdetect.getVirus(checkForAll[i].getKey(), checkForAll[i].getID()) == checkForAll[i]);
        } 
    }

    delete virus;
    return result;
}