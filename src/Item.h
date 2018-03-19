
//		Item.h

//		The Item class is an abstract base class for any
//		object which neesds to stored in a linked list.
//		The Itemlist class is a controller which gives most
//		of the functionality of the linked list. This implementation 
//		allows links to be stored with the Item instead of separately.
//
//		Nov. 6, 2001 - Added a double link and modified functions

#ifndef ITEM_H
#define ITEM_H
#include <stdlib.h>

class Item
{
	protected:
		int	n;				// Name (identification number) of object
		int index;			// Index position of Item (position in list)
		Item* nextItem;		// Pointer to next item in list
		Item* prevItem;		// Pointer to previous item in list
		
	public:
		Item(int num = 1, Item *next = NULL, Item *prev = NULL, int ind = -1)
		{ n=num; nextItem = next; index = ind;}
		int getn() const {return n;}
		void setn(int newn) {n=newn;}
		int getIndex() const {return index;}
		void setIndex(int newi) {index=newi;}
		Item* getNextOne() {return nextItem;}
		void setNextOne(Item* nextOne) {nextItem = nextOne;}
		Item* getPrevOne() {return prevItem;}
		void setPrevOne(Item* prevOne) {prevItem = prevOne;}
};

#endif
