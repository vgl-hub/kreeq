#ifndef KTREE_H
#define KTREE_H

class Knode {
    
    unsigned short int height;
    char* letter = NULL;
    Knode *A = NULL, *C = NULL, *G = NULL, *T = NULL;

public:
    
    Knode(unsigned short int height, char* letter) : height(height), letter(letter) {};
    
    void link(Knode* ptr);
    
    Knode* contains(char c);
    
    friend class Ktree;
    
};

class Ktree {
    
    Knode knodeRoot;
    unsigned short int KtreeH = 0;
    
public:
    
    Ktree (std::string* str, unsigned short int k);
    
    void print2D(Knode* parent, int space);
    
    void printKtree(Knode* root);
    
    void addChild(Knode* parent, unsigned long long int pos, unsigned short int height, char* c);
    
    void addKmer(char* c);
    
};

#endif /* KTREE_H */


