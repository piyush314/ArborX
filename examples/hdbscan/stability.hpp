#pragma once 


#include <map>
#include <vector>
#include <iostream>

typedef std::pair<int, int> score_t;

class stability_t
{
    public: 
    std::map<int, int> svect; 
    // default constructor 
    stability_t()
    {
        return ; 
    }

    // iterators to the map
    std::map<int, int>::iterator begin()
    {
        return svect.begin(); 
    }
    std::map<int, int>::iterator end()
    {
        return svect.end(); 
    }
    
    void sUpdate(int key, int val)
    {
        if(svect.find(key) == svect.end())
        {
            svect[key] = val; 
        }
        else
        {
            svect[key] += val; 
        }
        
    }
    // + operator overload with int
    stability_t operator+(int a)
    {
        stability_t temp = *this; 
        temp.sUpdate(a, 1);
         
        return temp; 
    }

    // + operator overload with stability_t
    stability_t operator+(stability_t& a)
    {
        stability_t temp = *this; 
        for(auto it = a.svect.begin(); it != a.svect.end(); it++)
        {
            temp.sUpdate(it->first, it->second);
            // temp.svect[it->first] += it->second;
        }
        return temp; 
    }

    // + operator overload with pair 
    stability_t operator+(std::pair<int, int> a)
    {
        stability_t temp = *this; 
        temp.sUpdate(a.first, a.second);
         
        return temp; 
    }

    // - operator overload with int
    stability_t operator-(int a)
    {
        stability_t temp = *this; 
        temp.sUpdate(a, -1);
        
        return temp; 
    }

    // - operator overload with stability_t
    stability_t operator-(stability_t& a)
    {
        stability_t temp = *this; 
        for(auto it = a.svect.begin(); it != a.svect.end(); it++)
        {
            temp.sUpdate(it->first, -it->second);
            
        }
        return temp; 
    }

    // - operator overload with pair
    stability_t operator-(std::pair<int, int> a)
    {
        stability_t temp = *this; 
        temp.sUpdate(a.first, -a.second);
        
        return temp; 
    }

    void print()
    {
        for(auto it = svect.begin(); it != svect.end(); it++)
        {
            if(it->second != 0)
            {
                std::cout << it->first << " " << it->second << "; ";
            }
            // std::cout << it->first << " " << it->second << std::endl;
        }
    }

    bool isZero()
    {
        for(auto it = svect.begin(); it != svect.end(); it++)
        {
            if(it->second != 0)
            {
                return false; 
            }
        }
        return true; 
    }

    void setZero()
    {
        for(auto it = svect.begin(); it != svect.end(); it++)
        {
            it->second = 0; 
        }
         
    }

    int getScore()
    {
        int score = 0; 
        for(auto it = svect.begin(); it != svect.end(); it++)
        {
            score += it->second; 
        }
        return score; 
    }

};


