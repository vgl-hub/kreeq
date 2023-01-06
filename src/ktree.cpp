#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <math.h>

#include <parallel_hashmap/phmap.h>

#include "bed.h"
#include "struct.h"
#include "global.h"
#include "uid-generator.h"
#include "gfa-lines.h"
#include "gfa.h"

#include "ktree.h"

void Ktree::DFS(Knode* current) {
    
    if(current == NULL)
        return;
    
    if (current->children[0] == NULL &&
        current->children[1] == NULL &&
        current->children[2] == NULL &&
        current->children[3] == NULL) {
        ++totKmersUnique;
        return;
        
    }
    
    // Process right children first
    DFS(current->children[0]);
    DFS(current->children[1]);
    
    // Process left children
    DFS(current->children[2]);
    DFS(current->children[3]);
    
    return;
    
}


void Ktree::print2D(Knode* current, int space) {
    
    if (current == NULL)
        return;
    
    space += 3;
    
    // Process right children first
    print2D(current->children[0], space);
    print2D(current->children[1], space);
    
    // print extended tree
    unsigned int height = current->height == 0 ? 1 : current->height;
    printf("\n%*s%.*s\n", space-3, "", height, current->letter);
    // print one level
//    printf("\n%*s%c\n", space-3, "", *(current->letter));
    
    // Process left children
    print2D(current->children[2], space);
    print2D(current->children[3], space);
    
    return;
    
}

void Ktree::printKtree(Knode* root) {

    print2D(root, 0);
    
}

//void Ktree::addChild(Knode* current, unsigned long long int pos, unsigned short int height, char* c) {
//
//    if (height>=KtreeH)
//        return;
//
//    switch (*(c+pos)) {
//
//        case 'A':
//
//            if (current->children[0] == NULL)
//                current-> = new Knode(c+pos);
//            addChild(current->A, ++pos, ++height, c);
//            break;
//
//        case 'C':
//
//            if (current->C == NULL)
//                current->C = new Knode(c+pos);
//            addChild(current->C, ++pos, ++height, c);
//            break;
//
//        case 'G':
//
//            if (current->G == NULL)
//                current->G = new Knode(c+pos);
//            addChild(current->G, ++pos, ++height, c);
//            break;
//
//        case 'T':
//
//            if (current->T == NULL)
//                current->T = new Knode(c+pos);
//            addChild(current->T, ++pos, ++height, c);
//
//    }
//
//    return;
//
//}

//void Ktree::addKmer(char* c) {
//    addChild(knodeRoot, 0, 0, c);
//}

void Ktree::addKmer(unsigned char* c) {
    
    Knode *parent = knodeRoot, *current = NULL;
    
    unsigned short int height = 0; // height is the offset in the node, pos is the height we are at
        
    current = parent->children[ctoi[*(c+height)]];
    
    for(unsigned short int pos = 0; pos < KtreeH; pos++) {
        
        if (pos+1 == KtreeH)
            return;
        
        if(current == NULL || *(current->letter+height++) != *(c++))
            break;
        
//        std::cout<<*(current->letter+height)<<"="<<*(c)<<" height is: "<<height<<" pos is: "<<pos<<"current height is: "<<current->height<<std::endl;

        
//        std::cout<<*(current->letter+height)<<"="<<*(c)<<" height is: "<<height<<" pos is: "<<pos<<"current height is: "<<current->height<<std::endl;
        
        if(current->height == height) {
            
            parent = current;
            
            current = current->children[ctoi[*(c)]];
            
            height = 0;
            
        }
        
    }
    
    if(current == NULL) {
        
        parent->children[ctoi[*(c)]] = new Knode(c);
        current = parent->children[ctoi[*(c)]];
        
    }else{
        
//        std::cout<<*(current->letter+height)<<"!="<<*(c)<<" height is: "<<height<<std::endl;
        
        if(current->height > height) { // case we need to branch and inherit children
            
            Knode* child = new Knode(current->letter+height); // create a new child node for current, the child starts at the branch
            
            memcpy(child->children, current->children, sizeof(current->children)); // the new child inherits the parent's children
            
            for (size_t i = 0; i < 4; i++) {current->children[i] = NULL;}// erase current childrens
            
            current->children[ctoi[*(current->letter+height)]] = child; // make the new child a children of current
            
//            std::cout<<"created branch A: "<<*current->children[ctoi[*(current->letter+height)]]->letter<<std::endl;
            
            if (current->children[ctoi[*(c)]] == NULL)
                current->children[ctoi[*(c)]] = new Knode(c); // new node

//            std::cout<<"created branch B: "<<*current->children[ctoi[*(c)]]->letter<<std::endl;
            
        }else{
            
            if (current->children[ctoi[*(current->letter+height)]] == NULL)
                current->children[ctoi[*(current->letter+height)]] = new Knode(current->letter+height); // original node
            
//            std::cout<<"created branch C: "<<*current->children[ctoi[*(current->letter+height)]]->letter<<std::endl;
            
            if (current->children[ctoi[*(c)]] == NULL)
                current->children[ctoi[*(c)]] = new Knode(c); // new node

//            std::cout<<"created branch D: "<<*current->children[ctoi[*(c)]]->letter<<std::endl;
            
        }
        
        current->height = height;
        
    }
    
    ++totKmersUnique;
    
}

Ktree::Ktree(InSequences& inSequences, unsigned short int k) {
    
    KtreeH = k;
    
    lg.verbose("Started ktree construction");
    
    knodeRoot = new Knode(new unsigned char('0'));
    
    std::vector<InSegment*>* segments = inSequences.getInSegments();
    
    for (InSegment* segment : *segments) {
        
        if (segment->getSegmentLen()<k) {
            lg.verbose("Segment " + segment->getSeqHeader() + " shorted thank k. skipping");
            continue;
        }
        
        unsigned long long int len = segment->getSegmentLen()-k+1;

        unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
        
        for (unsigned long long int c = 0; c<len; ++c) {

//            printf("adding kmer: %.*s\n", k, first+c);
            addKmer(first+c);
//            std::cout<<"Unique kmers: "<<totKmersUnique<<std::endl;
            
//            printKtree(knodeRoot);
            
        }
        
        totKmers += len;
        
        lg.verbose("Processed segment: " + segment->getSeqHeader());
        
    }
            
    printKtree(knodeRoot);
    
    std::cout<<"Total kmers: "<<totKmers<<std::endl;
    std::cout<<"Unique kmers: "<<totKmersUnique<<std::endl;
    
}

void Ktree::delKnodeRecurse(Knode* current) {
    
    if (current == NULL)
        return;

    delKnodeRecurse(current->children[3]);
    delKnodeRecurse(current->children[2]);
    delKnodeRecurse(current->children[1]);
    delKnodeRecurse(current->children[0]);
    
    delete current;
    
    return;
    
}

Ktree::~Ktree() {
    
    delete knodeRoot->letter;
//    delKnodeRecurse(knodeRoot);
    
}
