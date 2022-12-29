#ifndef KTREE_H
#define KTREE_H

class Knode {
    
    unsigned short int height;
    char* letter = NULL;
    Knode *A = NULL, *C = NULL, *G = NULL, *T = NULL;

public:
    
    Knode(unsigned short int height, char* letter) : height(height), letter(letter) {};
    
    Knode* contains(char c);
    
    friend class Ktree;
    
};

class Ktree {
    
    Knode* knodeRoot = NULL;
    unsigned short int KtreeH = 0;
    
public:
    
    Ktree(InSequences& inSequences, unsigned short int k);
    
    ~Ktree();
    
    void delKnodeRecurse(Knode* current);
    
    void print2D(Knode* current, int space);
    
    void printKtree(Knode* root);
    
    void addChild(Knode* current, unsigned long long int pos, unsigned short int height, char* c);
    
    void addKmer(char* c);
    
};

#endif /* KTREE_H */


