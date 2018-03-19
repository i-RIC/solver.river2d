#ifndef FILELOCKER_IMPL_H
#define FILELOCKER_IMPL_H

#include "../filelocker.h"

#ifdef _MSC_VER
	#include <windows.h>
#endif

class FileLocker::Impl
{
public:
	Impl(const std::string& filename);

#ifdef _MSC_VER
	HANDLE m_handle;
#else
	int m_fd;
#endif
	std::string m_filename;
	bool m_isLocked;
};

#endif // FILELOCKER_IMPL_H
