// utility.h
// Consider using boost library instead of using this header file.

#ifndef utility_h
#define utility_h

#include <string>

const std::string whiteSpaces( " \f\n\r\t\v" );

// ========================================
// Utility realted to string 
// Might be replaced by using boost/algorithm/string/trim.hpp

void trim_left( std::string& str,
      const std::string& trimChars = whiteSpaces )
{
   std::string::size_type pos = str.find_first_not_of( trimChars );
   str.erase( 0, pos );
}

void trim_right( std::string& str,
      const std::string& trimChars = whiteSpaces )
{
   std::string::size_type pos = str.find_last_not_of( trimChars );
   str.erase( pos + 1 );    
}

void trim( std::string& str, const std::string& trimChars = whiteSpaces )
{
   trim_left( str, trimChars );
   trim_right( str, trimChars );
}

// Distribute jobs through processors
// divide domain offset <= x < totalLength + offset
// Note that if mpi_size > totalLength && mpi_rank >= totalLength, 
// then istart = totalLength, iend = totalLength
void ParaRange(int totalLength, int offset, int mpi_size, int mpi_rank, int *istart, int *iend)
{
    int iwork1, iwork2;

    iwork1 = totalLength / mpi_size;
    iwork2 = totalLength % mpi_size;

    *istart = iwork1 * mpi_rank + min(iwork2, mpi_rank) + offset;
    *iend = *istart + iwork1;
    if (iwork2 > mpi_rank)
        *iend = *iend + 1;
}

#endif