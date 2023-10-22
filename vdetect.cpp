/***************************
 ** File:    vdetect.cpp
 ** Project: CMSC 341, proj4, Spring 2023
 ** Author:  William Whatley
 ** Date:    05/04/2023
 ** Section: Marron
 ** Email:   w59@umbc.edu
 **
 ** This file defines the functions laid out in vdetect.h. Every function is initialized here
 **
 ** All functions were created by William Whatley excluding isPrime, getNextPrime, and dump, whch was provided in the project
*****************************/

// CMSC 341 - Spring 2023 - Project 4
#include "vdetect.h"

//overloaded constructor
//constructs an empty object
VDetect::VDetect(int size, hash_fn hash, prob_t probing = DEFPOLCY){
    //if statement checks to see what the size integer sent is
    //if below MINPRIME, set equal to MINPRIME
    if (size <= MINPRIME){
        size = MINPRIME;
    }
    //else if size is greater than or equal to maxprime, set to maxprime
    else if (size >= MAXPRIME){
        size = MAXPRIME;
    }
    //else, checks to see if prime number
    else{
        //if it is not a prime number, findNextPrime is called so that the next prime number is found
        if (!isPrime(size)){
            size = findNextPrime(size);
        }
    }

    //other member variables are initialized
    m_currentSize = m_oldSize = 0; //sizes set to zero
    m_currNumDeleted = m_oldNumDeleted = 0; //numbers deleted set to zero
    m_currentCap = size; //currentCap is set to the size given
    m_oldCap = 0;
    m_oldProbing = NONE;
    m_currentTable = new Virus[size]; //currentTable is a new virus of size size 
    m_hash = hash;
    m_newPolicy = m_currProbing = probing; //newPolicy and currProbing set to probing
    m_oldTable = nullptr;
}

//destructor
//Deletes the old table and new table, and sets them to nullptr
VDetect::~VDetect(){
    delete []m_currentTable;
    delete []m_oldTable;
    m_currentTable = nullptr;
    m_oldTable = nullptr;

    //everything is reset to zero or default values
    m_currentSize = m_oldSize = 0; //sizes set to zero
    m_currNumDeleted = m_oldNumDeleted = 0; //numbers deleted set to zero
    m_currentCap = m_oldCap = 0;
    m_sizeAtStart = 0;
    m_newPolicy = m_currProbing = m_oldProbing = NONE;
}

//changeProbPolicy
//changes the probing policy for the hash table
void VDetect::changeProbPolicy(prob_t policy){
    //if statement checks to ensure that a valid policy is sent, and returns if it isn't
    if ((policy != NONE) && (policy != QUADRATIC) && (policy != DOUBLEHASH)){
        return;
    }
    m_newPolicy = policy; //policy is assigned to m_newPolicy
}

//insert
//inserts a virus and tries to find a hash for it
bool VDetect::insert(Virus virus){
    //if statement checks to ensure that the id is valid and that the table isn't at capacity AND also checks to ensure the key is valid
    if (((virus.getID() >= MINID) && (virus.getID() <= MAXID)) && (m_currentSize < (m_currentCap)) && (int(virus.getKey().length()) <= 5)){
        //checkVirus set to virus with the same key and ID
        Virus* checkVirus = new Virus(virus.getKey(), virus.getID());

        //getVirus used by sending the key of the current virus to see if it exists within the list
        *checkVirus = getVirus(virus.m_key, virus.m_id);

        //if statement checks to ensure that the virus is not already within the hash table by seeing if checkVirus is empty
        if (*checkVirus == EMPTY){
            int hash = m_hash(virus.getKey()); //the virus is hashed

            int hashNum = (hash % m_currentCap); //indexNum is the hash function divided by the current cap

            int indexNum = hashNum; //indexNum declared and initialized

            //if statement checks if the current spot is not empty and cannot be taken over, then determines what to do based on policy
            if (!(m_currentTable[indexNum] == EMPTY) && !(m_currentTable[indexNum] == DELETED)){
                //if statement checks which probing method is used
                if (m_currProbing == NONE){
                    //if the currentTable's position is occupied, which it means that the currentTable is neither EMPTY nor DELETED, that means it is occupied
                    delete checkVirus;
                    return false;
                }
                //else, checks if currProbing is a double hash
                else if(m_currProbing == DOUBLEHASH){
                    int i = 0; //i declared and initialized

                    //for loop goes through until m_capacity is reached to see if there is an open spot in the hashtable
                    for (; i < m_currentCap; i++){
                        indexNum = (hashNum + i * (11-(m_hash(virus.getKey()) % 11))) % m_currentCap;
                        //if a spot is empty or deleted, it means it is free
                        if ((m_currentTable[indexNum] == EMPTY) || (m_currentTable[indexNum] == DELETED)){
                            break; //for loop breaks early to indciate that a spot is found
                        }
                    }
                    //if i reaches current capacity, it means that it is unable to be hashed, thus false is returned
                    if (i == m_currentCap){
                        delete checkVirus;
                        return false;
                    }
                }
                //else, it is quadratic probing
                else{
                    int i = 0;
                    //for loop goes until m_currentCap is reached to see if an open probing spot exists
                    for (; i < (m_currentCap/2); i++){
                        indexNum = (hashNum + (i * i)) % m_currentCap;
                        //if a spot is empty or deleted, means open spot exists
                        if ((m_currentTable[indexNum] == EMPTY) || (m_currentTable[indexNum] == DELETED)){
                            break; //for loop breaks early to indciate that a spot is found
                        }
                    }

                    //if it reaches the end of the for loop and i = m_currentCap/2, it means that there's no spots available
                    if (i == (m_currentCap/2)){
                        delete checkVirus;
                        return false;
                    }
                }
            }

            //if all passes up to this point, that means m_currentTable spot can be used for checkVirus, thus it is saved
            m_currentTable[indexNum] = virus;

            ++m_currentSize; //size increased by one
            //if statement checks to see what the load factor is, and if it's greater than 0.5, the table is rehashed
            if (lambda() > 0.5){
                rehash();
                m_sizeAtStart = (m_oldSize - m_oldNumDeleted); //size at start is oldSize - oldNuMDeleted
            }

            transferData(); //transfer data called to transfer data over

            delete checkVirus; //checkVirus deletedW
            return true; //true returned
        }
        //else, if it's found not empty, checkVirus is deleted and false is returned
        else{
            delete checkVirus;
            return false;
        }
    }
    else{
        return false;
    }
}

//removes
//removes a Virus from the hash
bool VDetect::remove(Virus virus){
    //if statement checks to ensure that the ID is valid and that the string length of the virus sent is valid
    if (((virus.getID() >= MINID) && (virus.getID() <= MAXID)) && (int (virus.getKey().length()) <= 5)){
        Virus* checkVirus = new Virus(virus.getKey(), virus.getID());

        //getVirus used by sending the key of the current virus to see if it exists within the list
        *checkVirus = getVirus(virus.getKey(), virus.getID());

        //if statement checks to ensure that checkVirus
        if (!(*checkVirus == EMPTY)){
            int hashNum = m_hash(virus.getKey()) % m_currentCap; //key is hashed and divided by currentCap to see where located
            int index = hashNum; //index is set to hashNum
            int i = 0;
            bool found = false;

            //first thing checked is the currentTable, where the virus is seen if it can be removed
            if (m_currentSize != 0){
                //while loop goes through and checks to see if the virus is within the hashtable
                while (i < m_currentCap){
                    //if statement checks which probing method is used
                    if (m_currProbing == NONE){
                    //it checks to see if the initial hash is checkVirus, and found is true or false depending on it
                        if (m_currentTable[index] == *checkVirus){
                            found = true; //found is true since it was found
                            break; //breaks loop
                        }
                        else{
                            break; //breaks loop
                        }
                    }
                    //else, checks if currProbing is a double hash
                    else if(m_currProbing == DOUBLEHASH){
                        //if i != 0, then it changes the index
                        if (i != 0){
                            index = (hashNum + (i * (11-(m_hash(virus.getKey()) % 11)))) % m_currentCap;

                            //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                            if (m_currentTable[index] == *checkVirus){
                                found = true;
                                break;
                            }

                            //if it is empty, that means there is nothing afterwards, thus we can break
                            if (m_currentTable[index] == EMPTY){
                                break;
                            }
                        }
                        else{
                            //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                            if (m_currentTable[index] == *checkVirus){
                                found = true;
                                break;
                            }
                        }
                    }
                    //else, it is quadratic probing
                    else{
                        //if i != 0, then it changes the index
                        if (i != 0){

                            //of i is m_currentCap/2, then it breaks since doublehash has a limit of m_currentCap/2
                            if (i == (m_currentCap/2)){
                                break;
                            }
                            index = (hashNum + (i * i)) % m_currentCap;

                            //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                            if (m_currentTable[index] == *checkVirus){
                                found = true;
                                break;
                            }

                            //if it is empty, that means there is nothing afterwards, thus we can break
                            if (m_currentTable[index] == EMPTY){
                                break;
                            }
                        }
                        //else, initial hash is checked
                        else{
                            //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                            if (m_currentTable[index] == *checkVirus){
                                found = true;
                                break;
                            }
                        }
                    }
                    ++i; //i size increased by one
                }
            }  


            //if it is found in the current table, it is deleted, and deleted ratio is considered
            if (found){
                m_currentTable[index] = DELETED; //it becomes deleted
                ++m_currNumDeleted; //the currentNumDeleted is increased by one
                
                //if the deleted ratio is greater than 0.8, or numDeleted that means the table needs to be rehashed
                if (deletedRatio() > 0.8){
                    rehash();
                    m_sizeAtStart = (m_oldSize - m_oldNumDeleted); //size at start is oldSize - oldNuMDeleted
                }

                transferData(); //transferData called to transfer data over

                delete checkVirus; //checkVirus deleted;

                return true; //true is returned
            }

            //else, if it's not in the newTable, it must be in the oldTable, so it is then searched depending on that

            hashNum = m_hash(virus.getKey()) % m_oldCap; //key is hashed and divided by oldCap to see where located
            index = hashNum; //index set to the new hash number
            i = 0;

            //while loop goes through the previous table and checks to see if it exists
            while (i < m_oldCap){
                //if statement checks which probing method is used
                if (m_oldProbing == NONE){
                    break; //if the old probing is none, that means that it must be the first one to show up on the table
                }
                //else, checks if currProbing is a double hash
                else if(m_oldProbing == DOUBLEHASH){
                    //if i != 0, then it changes the index
                    if (i != 0){
                        index = (hashNum + i * (11-(m_hash(virus.getKey()) % 11))) % m_oldCap;

                        //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                        if (m_oldTable[index] == *checkVirus){
                            break;
                        }
                    }
                    else{
                        //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                        if (m_oldTable[index] == *checkVirus){
                            break;
                        }
                    }
                }
                //else, it is quadratic probing
                else{
                    //if i != 0, then it changes the index
                    if (i != 0){
                        index = (hashNum + (i * i)) % m_oldCap;

                        //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                        if (m_oldTable[index] == *checkVirus){
                            break;
                        }
                    }
                    //else, initial hash is checked
                    else{
                        //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                        if (m_oldTable[index] == *checkVirus){
                            break;
                        }
                    }
                }
                ++i; //i size increased by one
            }

            m_oldTable[index] = DELETED; //it becomes deleted from the old table
            ++m_oldNumDeleted; //the oldNumDeleted is increased by one

            transferData(); //transferData is called if it is needed

            delete checkVirus;
            return true;
        }
        //else, if an empty is returned, that means it is not within the tables, thus false is returned
        else{
            delete checkVirus;
            return false;
        }
    }
    else{
        return false; //else, false is returned
    }
}

//getVirus
//returns the key and index of the virus function
Virus VDetect::getVirus(string key, int id) const{
    //if statement checks to ensure the key and ID sent are valid
    if (((id >= MINID) && (id <= MAXID)) && (int(key.length()) <= 5)){
        Virus* checkVirus = new Virus(key, id); //checkVirus is a pointer with the new virus data
        int hashNum = m_hash(key) % m_currentCap; //key is hashed and divided by currentCap to see where located
        int index = hashNum; //index is set to hashNum
        int i = 0;
        bool found = false;

        //if statement checks to ensure the current size is not zero
        if (m_currentSize != 0){
            //while loop goes through and checks to see if the virus is within the hashtable
            while (i < m_currentCap){
                //if statement checks which probing method is used
                if (m_currProbing == NONE){
                //it checks to see if the initial hash is checkVirus, and found is true or false depending on it
                    if (m_currentTable[index] == *checkVirus){
                        found = true; //found is true since it was found
                        break; //breaks loop
                    }
                    else{
                        break; //breaks loop
                    }
                }
                //else, checks if currProbing is a double hash
                else if(m_currProbing == DOUBLEHASH){
                    //if i != 0, then it changes the index
                    if (i != 0){
                        index = (hashNum + (i * (11-(m_hash(key) % 11)))) % m_currentCap;

                        //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                        if (m_currentTable[index] == *checkVirus){
                            found = true;
                            break;
                        }

                        //if it is empty, that means there is nothing afterwards, thus we can break
                        if (m_currentTable[index] == EMPTY){
                            break;
                        }
                    }
                    else{
                        //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                        if (m_currentTable[index] == *checkVirus){
                            found = true;
                            break;
                        }
                    }
                }
                //else, it is quadratic probing
                else{
                    //if i != 0, then it changes the index
                    if (i != 0){
                        index = (hashNum + (i * i)) % m_currentCap;

                        //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                        if (m_currentTable[index] == *checkVirus){
                            found = true;
                            break;
                        }

                        //else if there is an empty index, that means that there is nothing afterwards, thus we can break
                        if (m_currentTable[index] == EMPTY){
                            break;
                        }

                        //if we hit m_currentCap/2 - 1, that means there is nothing more we can find, thus we break
                        if (i == (m_currentCap/2 - 1)){
                            break;
                        }
                    }
                    //else, initial hash is checked
                    else{
                        //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                        if (m_currentTable[index] == *checkVirus){
                            found = true;
                            break;
                        }
                    }
                }
                ++i; //i size increased by one
            }
        }
        
        //if statement continues if the virus is not found still and also checks that the oldTable isn't nullptr
        if (!(found) && (m_oldTable != nullptr)){
            //if statement checks to ensure the size of old table isn't zero
            if (m_oldSize != 0){
                hashNum = m_hash(key) % m_oldCap; //key is hashed and divided by oldCap to see where located
                index = hashNum; //index set to the new hash number
                i = 0;

                //while loop goes through the previous table and checks to see if it exists
                while (i < m_oldCap){
                    //if statement checks which probing method is used
                    if (m_oldProbing == NONE){
                    //it checks to see if the initial hash is checkVirus, and found is true or false depending on it
                        if (m_oldTable[index] == *checkVirus){
                            found = true;
                            break;
                        }
                        else{
                            break; //loop breaks
                        }
                    }
                    //else, checks if currProbing is a double hash
                    else if(m_oldProbing == DOUBLEHASH){
                        //if i != 0, then it changes the index
                        if (i != 0){
                            index = (hashNum + i * (11-(m_hash(key) % 11))) % m_oldCap;

                            //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                            if (m_oldTable[index] == *checkVirus){
                                found = true;
                                break;
                            }

                            //if we hit empty, we break since there will be nothing else after it
                            if (m_oldTable[index] == EMPTY){
                                break;
                            }
                        }
                        else{
                            //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                            if (m_oldTable[index] == *checkVirus){
                                found = true;
                                break;
                            }
                        }
                    }
                    //else, it is quadratic probing
                    else{
                        //if i != 0, then it changes the index
                        if (i != 0){
                            //if i is  m_oldCap/2, then it breaks since doublehash has that as a limit
                            if (i == (m_oldCap/2)){
                                break;
                            }
                            index = (hashNum + (i * i)) % m_oldCap;

                            //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                            if (m_oldTable[index] == *checkVirus){
                                found = true;
                                break;
                            }
                                                        
                                                        
                            //if we hit empty, we break since there will be nothing else after it
                            if (m_oldTable[index] == EMPTY){
                                break;
                            }

                        }
                        //else, initial hash is checked
                        else{
                            //if statement checks that it is occupied by the correct virus, and then breaks the loop if it is
                            if (m_oldTable[index] == *checkVirus){
                                found = true;
                                break;
                            }
                        }
                    }
                ++i; //i size increased by one
                }
            }
        }

        //if statement checks if it is found
        if (found){
            //virus is declared and assigned to virus check
            Virus virus(key, id);
            virus = *checkVirus;
            delete checkVirus; //virus check is deleted
            return virus; //virus is returned
        }
        else{
            delete checkVirus;
            return EMPTY; //empty is returned to indicate that it isn't found
        }
    }
    else{
        return EMPTY; //else, empty is returned
    }
}

//lambda
//calculates the load factor of the function
float VDetect::lambda() const {
    return float(m_currentSize)/float(m_currentCap);
}

//deletedRatio
//returns curr number deleted by the current size
float VDetect::deletedRatio() const {
    return float(m_currNumDeleted)/float(m_currentSize);
}

void VDetect::dump() const {
    cout << "Dump for the current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for the old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

bool VDetect::isPrime(int number){
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int VDetect::findNextPrime(int current){
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) { 
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0) 
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}

ostream& operator<<(ostream& sout, const Virus &virus ) {
    if (!virus.m_key.empty())
        sout << virus.m_key << " (ID " << virus.m_id << ")";
    else
        sout << "";
  return sout;
}

bool operator==(const Virus& lhs, const Virus& rhs){
    return ((lhs.m_key == rhs.m_key) && (lhs.m_id == rhs.m_id));
}

void VDetect::rehash(){
    //if oldTable isn't a nullptr, it means it already exists, thus nothing can be done
    if (m_oldTable != nullptr){
        return;
    }

    //first, all member variables are swapped from the old table to the new table
    m_oldTable = m_currentTable;
    m_oldNumDeleted = m_currNumDeleted;
    m_oldSize = m_currentSize;
    m_oldCap = m_currentCap;
    m_oldProbing = m_currProbing;

    //next, the member variables of the currentTable are changed
    m_currentSize = 0;
    m_currNumDeleted = 0;
    m_currProbing = m_newPolicy;

    //m_currentCapacity has a new size created
    int newCap = 4 * (m_oldSize - m_oldNumDeleted);


    //newCap set to findNextPrime, then that is checked if it is valid
    newCap = findNextPrime(newCap);
    //if it's less than MINPRIME, new capacity is minprime
    if (newCap < MINPRIME){
        m_currentCap = MINPRIME;
    }
    //else, if it's greater than MAXPRIME, it becomes maxprime
    else if (newCap > MAXPRIME){
        m_currentCap = MAXPRIME;
    }
    //else, it is set to the newCap if it is a good prime number
    else{
        m_currentCap = newCap;
    }

    //the new table is created with m_currentCap as the size
    m_currentTable = new Virus[m_currentCap];
}

//transferData
//transfers data from the old table to the new table if it's required
void VDetect::transferData(){
    //if m_oldTable is a nullptr, means that no transferring is required, alongside if the current table is at capacity
    if ((m_oldTable == nullptr) || (m_currentSize == m_currentCap)){
        return;
    }

    //iterary operator does the starting point  
    int endBound = (m_sizeAtStart / 4);
    int count = 0;

    //if the endBound is 0, that means there's between 1-3 data points still on the table, so it is reset
    if (endBound == 0){
        endBound = (m_oldSize - m_oldNumDeleted);
    }
    
    //for loop goes through until it hits the starting point and adds 25% of the data
    for (int i = 0; i < m_oldCap; i++){
        //if the spot isn't empty
        if (!(m_oldTable[i] == EMPTY) && !(m_oldTable[i] == DELETED)){
            //virus ptr set to oldTable[i], then oldTable[i] set to deleted so that insertion can be tested
            if (uniqueInsert(m_oldTable[i])){
                m_oldTable[i] = DELETED;
                ++m_oldNumDeleted; //old num deleted increased by one
                ++count;
            }
        }

        if ((count == endBound) || ((m_oldSize - m_oldNumDeleted) == 0)){
            break;
        }
    }


    //if the old size becomes zero, everything is reset to default in the old table
    if ((m_oldSize - m_oldNumDeleted) <= 0){
        delete []m_oldTable;
        m_oldTable = nullptr;
        m_oldNumDeleted = m_oldSize = 0; 
        m_oldCap = 0;
        m_oldProbing = NONE;
        m_sizeAtStart = 0;
    }
}

bool VDetect::uniqueInsert(Virus virus){
    Virus* checkVirus = new Virus();

    *checkVirus = virus;

    int hash = m_hash(virus.getKey()); //the virus is hashed

    int hashNum = (hash % m_currentCap); //indexNum is the hash function divided by the current cap

    int indexNum = hashNum; //indexNum declared and initialized

    //if statement checks if the current spot is not empty and cannot be taken over, then determines what to do based on policy
    if (!(m_currentTable[indexNum] == EMPTY) && !(m_currentTable[indexNum] == DELETED)){
        //if statement checks which probing method is used
        if (m_currProbing == NONE){
            //if the currentTable's position is occupied, which it means that the currentTable is neither EMPTY nor DELETED, that means it is occupied
            delete checkVirus;
            return false;
        }
        //else, checks if currProbing is a double hash
        else if(m_currProbing == DOUBLEHASH){
            int i = 0; //i declared and initialized

            //for loop goes through until m_capacity is reached to see if there is an open spot in the hashtable
            for (; i < m_currentCap; i++){
                indexNum = (hashNum + i * (11-(m_hash(virus.getKey()) % 11))) % m_currentCap;
                //if a spot is empty or deleted, it means it is free
                if ((m_currentTable[indexNum] == EMPTY) || (m_currentTable[indexNum] == DELETED)){
                    break; //for loop breaks early to indciate that a spot is found
                }
            }
            //if i reaches current capacity, it means that it is unable to be hashed, thus false is returned
            if (i == m_currentCap){
                delete checkVirus;
                return false;
            }
        }
        //else, it is quadratic probing
        else{
            int i = 0;
            //for loop goes until m_currentCap is reached to see if an open probing spot exists
            for (; i < (m_currentCap/2); i++){
                indexNum = (hashNum + (i * i)) % m_currentCap;
                //if a spot is empty or deleted, means open spot exists
                if ((m_currentTable[indexNum] == EMPTY) || (m_currentTable[indexNum] == DELETED)){
                    break; //for loop breaks early to indciate that a spot is found
                }
            }

            //if it reaches the end of the for loop and i = m_currentCap/2, it means that there's no spots available
            if (i >= (m_currentCap/2)){
                delete checkVirus;
                return false;
            }
        }
    }

    //if all passes up to this point, that means m_currentTable spot can be used for checkVirus, thus it is saved
    m_currentTable[indexNum] = virus;

    ++m_currentSize; //size increased by one

    delete checkVirus; //checkVirus deletedW
    return true; //true returned
}