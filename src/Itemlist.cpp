//
//		ItemList.cpp
//

#include "Itemlist.h"
#include <iostream>
using namespace std;
	
ItemList::ItemList()
{
	numItems = 0;
	theFirstItem = NULL;
	theCurrentItem = NULL;
	itemIndex = NULL;
	indexState = notValid;
}

ItemList::~ItemList()
{
	delete [] itemIndex;
}

void ItemList::emptyList()
{
	firstItem();
	while(numItems>0)
		deleteCurrentItem();
	delete [] itemIndex;
	itemIndex = NULL;
	indexState = notValid;
}

Item* ItemList::i(int n)
{
	if(indexState == notValid) {
		buildIndex();
	}
	if((n <= numItems) && (n > 0))
		return (itemIndex[n-1]);
	else
		return NULL;
}

int ItemList::buildIndex()
{
	int index;

	delete [] itemIndex;
	if(numItems > 0) {
		itemIndex = new Item* [numItems];
		itemIndex[0] = firstItem();
		theFirstItem->setIndex(0);
		for(index=1;index<numItems;index++){
			itemIndex[index] = nextItem();
			theCurrentItem->setIndex(index);
		}
		indexState = valid;
	}
	else {
		itemIndex = NULL;
		indexState = notValid;
	}
	return numItems;
}

Item* ItemList::firstItem()
{
	if(numItems >0){
		theCurrentItem = theFirstItem;
		return(theFirstItem);
	}
	else
		return(NULL);
}							
			
Item* ItemList::currentItem()
{
	if(numItems >0){
		return(theCurrentItem);
	}
	else
		return(NULL);
}							
			
Item* ItemList::setCurrentItem(Item *itemP)
{
	if(itemP != NULL)
		theCurrentItem = itemP;
	else
		theCurrentItem = theFirstItem;
		
	return theCurrentItem;
}							
			
Item* ItemList::nextItem()
{
	if(numItems > 0) {
		if(theCurrentItem->getNextOne() != NULL){
			theCurrentItem = theCurrentItem->getNextOne();
			return (theCurrentItem);
		}
		else
			return(NULL);
	}
	else
		return(NULL);
}
			
Item* ItemList::n(int name)
{
	if(numItems > 0) {
		theCurrentItem = theFirstItem;
		while(theCurrentItem->getn() != name){
			if(theCurrentItem->getNextOne() != NULL){
				theCurrentItem = theCurrentItem->getNextOne();
			}
			else
				return(NULL);
		}
		return(theCurrentItem);
	}
	else
		return(NULL);
}

Item* ItemList::lastItem()
{
	if(numItems > 0) {
		while(theCurrentItem->getNextOne() != NULL) {
			theCurrentItem = theCurrentItem->getNextOne();
		}
	}
	return (theCurrentItem);
}
		
Item* ItemList::appendItem(Item *theNewItem)
{
	if(numItems > 0) {
		while(theCurrentItem->getNextOne() != NULL) {
			theCurrentItem = theCurrentItem->getNextOne();
		}
	}
	return (insertItem(theNewItem));
}
		
Item* ItemList::insertItem(Item *theNewItem)
{	
	if(numItems > 0) {
		Item* theNextItem = theCurrentItem->getNextOne();
		theNewItem->setNextOne(theNextItem);
		theNewItem->setPrevOne(theCurrentItem);
		if(theNextItem != NULL)
			theNextItem->setPrevOne(theNewItem);
		theCurrentItem->setNextOne(theNewItem);
		theCurrentItem = theNewItem;
	}
	else {
		theNewItem->setNextOne(NULL);
		theNewItem->setPrevOne(NULL);
		theFirstItem = theNewItem;
		theCurrentItem = theNewItem;
	}
	numItems += 1;
	indexState = notValid;
	
	return(theNewItem);
}

Item* ItemList::push(Item *theNewItem)
{	
	if(numItems > 0) {
		theNewItem->setNextOne(theFirstItem);
		theNewItem->setPrevOne(NULL);
		theFirstItem = theNewItem;
		theCurrentItem = theNewItem;
	}
	else {
		theNewItem->setNextOne(NULL);
		theNewItem->setPrevOne(NULL);
		theFirstItem = theNewItem;
		theCurrentItem = theNewItem;
	}
	numItems += 1;
	indexState = notValid;
	
	return(theNewItem);
}

Item* ItemList::pop()
{	
	Item* popItemP;
	
	if(numItems > 0) {
		popItemP = theFirstItem;
		theFirstItem = popItemP->getNextOne();
		if(theFirstItem != NULL)
			theFirstItem->setPrevOne(NULL);
		theCurrentItem = theFirstItem;
		numItems -= 1;
		indexState = notValid;
		
		return popItemP;
	}
	else 	
		return NULL;
}

void ItemList::clearList()
{		
	Item* nextItemP;
	
	theCurrentItem = theFirstItem;
	
	while(numItems > 0) {
		nextItemP = theCurrentItem->getNextOne();
//		try{
			delete theCurrentItem;
//		}
//		catch(...){
			// ignore failure
//		}
		theCurrentItem = nextItemP;
		numItems -= 1;
	}
	theFirstItem = NULL;
	delete [] itemIndex;
	itemIndex = NULL;
	indexState = notValid;
}

void ItemList::catList(ItemList *otherList)
{
	if(numItems > 0) {
		while(theCurrentItem->getNextOne() != NULL) {
			theCurrentItem = theCurrentItem->getNextOne();
		}
		Item *otherFirstItem = otherList->firstItem();
		theCurrentItem->setNextOne(otherFirstItem);
		if(otherFirstItem != NULL)
			otherList->firstItem()->setPrevOne(theCurrentItem);
		theCurrentItem = theFirstItem;
	}
	else {
		theFirstItem = otherList->firstItem();
		theCurrentItem = theFirstItem;
	}
	numItems += otherList->numberOfItems();
	indexState = notValid;
}

Item* ItemList::deleteCurrentItem()
{
	Item *prevItem;

	prevItem = theFirstItem;

	if(theFirstItem != NULL) {
		if(theCurrentItem != theFirstItem) {
			while(prevItem->getNextOne() != theCurrentItem) {
				prevItem = prevItem->getNextOne();
				if(prevItem == NULL){
					return theCurrentItem;
				}
			}
			prevItem->setNextOne(theCurrentItem->getNextOne());
			if(prevItem->getNextOne() != NULL)
				theCurrentItem = prevItem->getNextOne();
			else
				theCurrentItem = prevItem;
		}
		else{
			theFirstItem = theCurrentItem->getNextOne();
			theCurrentItem = theFirstItem;
		}
		numItems -= 1;
	}
	indexState = notValid;
	return (theCurrentItem);
}

Item* ItemList::deleteItem(Item *iP)
{
	Item *prevItem = iP->getPrevOne();
	Item *nextItem = iP->getNextOne();

	if(theFirstItem == NULL)
		return NULL;

	if(iP == theFirstItem) {
		theFirstItem = nextItem;
		if(nextItem != NULL)
			nextItem->setPrevOne(NULL);
		numItems -= 1;
		indexState = notValid;
		if(iP == theCurrentItem){
			theCurrentItem = theFirstItem;
		}
		return theCurrentItem;
	}

	if(iP == prevItem->getNextOne()) {
		prevItem->setNextOne(nextItem);
		if(nextItem != NULL)
			nextItem->setPrevOne(prevItem);
		numItems -= 1;
		indexState = notValid;
		if(iP == theCurrentItem)
				theCurrentItem = prevItem;
	}
	return theCurrentItem;
}

