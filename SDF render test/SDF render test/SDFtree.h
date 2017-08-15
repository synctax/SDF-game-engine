//
//  SDFtree.h
//  SDF renderer with shader
//
//  Created by Ben Wyngaard on 8/1/17.
//  Copyright © 2017 Ben Wyngaard. All rights reserved.
//
#include <glm/vec3.hpp> // glm::vec3
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtc/matrix_transform.hpp>
#include <algorithm>


#ifndef SDFtree_h
#define SDFtree_h
namespace SDF{
    struct sdfNode{
        int nodeType;
        int subType;
        float args[3];
        int extra;
        int flags;
        sdfNode* childNodes[2];
        int children[2];
        int addr;
    };
    
    sdfNode fSphere(float radius, int matID,bool isBounding){
        return {1,1,{radius,0,0},matID,(int)isBounding,{nullptr,nullptr}};
    }
    sdfNode fBox(glm::vec3 size, int matID,bool isBounding){
        return {1,2,{size.x,size.y,size.z},matID,(int)isBounding,{nullptr,nullptr}};
    }
    sdfNode fPlane(glm::vec3 normal, int matID,bool isBounding){
        return {1,3,{normal.x,normal.y,normal.z},matID,(int)isBounding,{nullptr,nullptr}};
    }
    sdfNode fCapsule(float radius, float length, int matID,bool isBounding){
        return {1,4,{radius,length,0},matID,(int)isBounding,{nullptr,nullptr}};
    }
    sdfNode fCylinder(float radius, float height,int matID,bool isBounding){
        return {1,5,{radius,height,0},matID,(int)isBounding,{nullptr,nullptr}};
    }
    sdfNode fCone(float radius, float height, int matID,bool isBounding){
        return {1,6,{radius,height,0},matID,(int)isBounding,{nullptr,nullptr}};
    }
    sdfNode fTorus(float smallR, float bigR, int matID,bool isBounding){
        return {1,7,{smallR,bigR,0},matID,(int)isBounding,{nullptr,nullptr}};
    }
    
    
    sdfNode fOpUnion(sdfNode* a, sdfNode* b){
        return {2,1,{0,0,0},0,0,{a,b}};
    }
    sdfNode fOpIntersect(sdfNode* a, sdfNode* b){
        return {2,2,{0,0,0},0,0,{a,b}};
    }
    
    sdfNode fOpDifference(sdfNode* a, sdfNode* b){
        return {2,3,{0,0,0},0,0,{a,b}};
    }
    
    sdfNode pTranslate(sdfNode* a, glm::vec3 amount){
        return {3,1,{amount.x,amount.y,amount.z},0,0,{a,nullptr}};
    }
    
    int* decomposeNode(sdfNode* node){
        int ints[4];
        ints[0] = (int)node->nodeType;
        ints[0] += ((int)node->subType)<<8;
        ints[0] += ((int)node->extra)<<16;
        ints[0] += ((int)node->flags)<<24;
        ints[1] = *(int*)&node->args[0];
        ints[2] = *(int*)&node->args[1];
        ints[3] = *(int*)&node->args[2];
        return ints;
    }
    
    class sdfTree{
    private:
        sdfNode* head;
        
        int findNextOpen(){
            int cur = 0;
            while(this->buffer[cur] != 0){
                cur += 5;
            }
            return cur;
        }
        
        int fillBuffer(sdfNode* n, int index){
            this->maxIndex = std::max(index, this->maxIndex);
            this->buffer[index] = 1;
            int* ints = decomposeNode(n);
            if (n->childNodes[0] != nullptr) ints[4] = (this->fillBuffer(n->childNodes[0],index*2) << 16);
            if (n->childNodes[1] != nullptr) ints[4] = (this->fillBuffer(n->childNodes[1], (index*2)+1));
            this->buffer[index] = ints[0]; this->buffer[index+1] = ints[1]; this->buffer[index+2] = ints[2];
            this->buffer[index+3] = ints[3]; this->buffer[index+4] = ints[4];
            n->addr = index;
            return index;
        }
    public:
        int buffer[50], maxIndex = 0;
        
        sdfTree(sdfNode* head){
            this->head = head;
            head->addr = 0;
            this->fillBuffer(head, 0);
        }
        
        void updateNode(sdfNode* node){
            if (node->addr){
                int index = node->addr;
                int* ints = decomposeNode(node);
                this->buffer[index] = ints[0]; this->buffer[index+1] = ints[1]; this->buffer[index+2] = ints[2];
                this->buffer[index+3] = ints[3];
            }
        }
    };
    
}



#endif /* SDFtree_h */
