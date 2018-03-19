#ifndef FILELOCKER_H
#define FILELOCKER_H

#include "iriclib_global.h"

#include <string>

class IRICLIBDLL FileLocker
{
public:
	FileLocker(const std::string& filename);
	~FileLocker();

	bool lockFileExists() const;

	bool isLocked() const;

	bool lock();
	void unlock();

private:
	class Impl;
	Impl* impl;
};

#ifdef _DEBUG
	#include "private/filelocker_impl.h"
#endif // _DEBUG

#endif // FILELOCKER_H
