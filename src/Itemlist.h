//   ItemList Class

//	The ItemList Class is intended to provide a fairly generic
//	facility for linked list storage and access of a collection of objects.
//	The objects must derive from the Item class to contain 
//	the necessary list variables. A cast is necessary when getting items,
//	but not when putting them.
//	Most operations defined in Item List
//  affect the organization of the entries only.
//	The list store does not manage memory,
//	either allocation or freeing, of the items themselves. Thus, in clearing,
//  first delete the item, then the ItemList entry. Deleting the ItemList itself 
//  will delete all of the information about the list, but none of the items.
//
//

#ifndef ITEMLIST_H
#define ITEMLIST_H

#include "Item.h"

enum status {valid, notValid};
	
class ItemList
{
	
	protected:
		Item	*theFirstItem;
		Item	*theCurrentItem;
		int		numItems;
		Item	**itemIndex;
		status	indexState;
		
			
	public:
		ItemList();								//	construcor
		~ItemList();							//	destructor
		void emptyList();						//	does not delete items
		int numberOfItems() {return numItems;}		
		Item* i(int i);							//	return pointer to ith item
		int buildIndex();
		Item* firstItem();						//	get first item in list
		Item* currentItem();					//	get current item in list
		Item* setCurrentItem(Item *ItemP);		//	set current item in list
		Item* nextItem();						//	get next item in list
		Item* n(int name);						//	get item with n = name
		Item* lastItem();
		Item* appendItem(Item *theNewItem);		//	add new item to end of list
		Item* insertItem(Item *theNewItem);	//	insert new item after current item
		Item* push(Item *theNewItem);			//	put new item at start of list
		Item* pop();							//	take (remove) item from start of list
		void clearList();						//	delete (really delete) all items in list
		void catList(ItemList *otherList);		//	concatenate another list to end of list
		Item* deleteCurrentItem();				//	remove current item, next becomes current
		Item* deleteItem(Item *iP);		// delete item pointed at, return currentItem
};

#endif
