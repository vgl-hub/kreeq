#include <stdlib.h>
#include <string>
#include <iostream>
#include <math.h>

#include "global.h"

#include "ktree.h"

void Knode::link(Knode* knodePtr){
    
    switch(*(knodePtr->letter)) {
        case 'A':
            this->A = knodePtr;
            break;
        case 'C':
            this->C = knodePtr;
            break;
        case 'G':
            this->G = knodePtr;
            break;
        case 'T':
            this->T = knodePtr;
            break;
    }
}

Knode* Knode::contains(char c){
    
    switch(c) {
        case 'A':
            if (this->A != NULL)
                return this->A;
            break;
        case 'C':
            if (this->C != NULL)
                return this->C;
            break;
        case 'G':
            if (this->G != NULL)
                return this->G;
            break;
        case 'T':
            if (this->T != NULL)
                return this->T;
            break;
    }
    
    return NULL;
    
}

void Ktree::print2D(Knode* parent, int space) {

    if (parent == NULL)
        return;
 
    space += 3;
 
    // Process right children first
    print2D(parent->T, space);
    print2D(parent->G, space);
 
    printf("\n%*s%c\n", space-3, "", *(parent->letter));
 
    // Process left children
    print2D(parent->C, space);
    print2D(parent->A, space);
}

void Ktree::printKtree(Knode* root) {
    print2D(root, 0);
}

void Ktree::addChild(Knode* parent, unsigned long long int pos, unsigned short int height, char* c) {
    
    if (height>=KtreeH)
        return;
    
    switch (*(c+pos)) {
            
        case 'T':
            
            if (parent->T == NULL)
                parent->T = new Knode(height, c+pos);
            addChild(parent->T, ++pos, ++height, c);
            break;
            
        case 'G':
            
            if (parent->G == NULL)
                parent->G = new Knode(height, c+pos);
            addChild(parent->G, ++pos, ++height, c);
            break;
            
        case 'C':
            
            if (parent->C == NULL)
                parent->C = new Knode(height, c+pos);
            addChild(parent->C, ++pos, ++height, c);
            break;
            
        case 'A':
            
            if (parent->A == NULL)
                parent->A = new Knode(height, c+pos);
            addChild(parent->A, ++pos, ++height, c);
            
    }
    
    return;
}

void Ktree::addKmer(char* c) {
    addChild(&knodeRoot, 0, 0, c);
}

Ktree::Ktree (std::string* str, unsigned short int k) :
knodeRoot(0, new char('0')) {
    
    KtreeH = k;
    
    lg.verbose("Started ktree construction");
    
    unsigned long long int len = str->size()-k+1;

    for (unsigned short int c = 0; c<len; ++c) {
        
        addKmer(&(*str)[c]);
        
    }
    
//    knodePtr->A = new Knode(0, 'A');
//    knodePtr->A->C = new Knode(0, 'C');
//    knodePtr->A->C->C = new Knode(0, 'C');
//
//    knodePtr->C = new Knode(0, 'C');
//    knodePtr->C->C = new Knode(0, 'C');
//    knodePtr->C->C->T = new Knode(0, 'T');
//
//    knodePtr->C->T = new Knode(0, 'T');
//    knodePtr->C->T->G = new Knode(0, 'G');
//
//    knodePtr->T = new Knode(0, 'T');
//    knodePtr->T->G = new Knode(0, 'G');
//    knodePtr->T->G->C = new Knode(0, 'C');
//
//    knodePtr->G = new Knode(0, 'G');
//    knodePtr->G->C = new Knode(0, 'C');
//    knodePtr->G->C->C = new Knode(0, 'C');
//
//    knodePtr->T->G->A = new Knode(0, 'A');

    
//    Knode* knodePtr = &knodeRoot;
//
//    unsigned long long int len = str.size()-k+1;
//
//    for (unsigned short int c = 0; c<len; ++c) {
//
//        Knode* knodeCurr = knodeRoot.contains(str[c]);
//
//        if (knodeCurr != NULL)
//            continue;
//
//        Knode knode(0, str[c]);
//
//
//    }
            
    printKtree(&knodeRoot);
    
}
