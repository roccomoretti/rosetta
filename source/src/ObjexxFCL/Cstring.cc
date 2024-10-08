// Cstring: C String Wrapper
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// ObjexxFCL Headers
#include <ObjexxFCL/Cstring.hh>

// C++ Headers
#include <cctype>
#include <iostream>


namespace ObjexxFCL {


// Constants
char const SPACE( ' ' );
char const TAB( '\t' );
char const NULLC( '\000' );


// Cstring


/// @brief Has Any Character of a Cstring?
bool
Cstring::has_any_of( Cstring const & s ) const
{
	size_type const s_len( s.length() );
	for ( size_type i = 0; i < std::strlen( str_ ); ++i ) {
		for ( size_type j = 0; j < s_len; ++j ) {
			if ( str_[ i ] == s.str_[ j ] ) return true;
		}
	}
	return false; // No matches
}


/// @brief Has Any Character of a cstring?
bool
Cstring::has_any_of( c_cstring const & s ) const
{
	size_type const s_len( std::strlen( s ) );
	for ( size_type i = 0; i < std::strlen( str_ ); ++i ) {
		for ( size_type j = 0; j < s_len; ++j ) {
			if ( str_[ i ] == s[ j ] ) return true;
		}
	}
	return false; // No matches
}


/// @brief Has Any Character of a std::string?
bool
Cstring::has_any_of( std::string const & s ) const
{
	size_type const s_len( s.length() );
	for ( size_type i = 0; i < std::strlen( str_ ); ++i ) {
		for ( size_type j = 0; j < s_len; ++j ) {
			if ( str_[ i ] == s[ j ] ) return true;
		}
	}
	return false; // No matches
}


/// @brief Has a Character?
bool
Cstring::has_any_of( char const c ) const
{
	for ( size_type i = 0; i < std::strlen( str_ ); ++i ) {
		if ( str_[ i ] == c ) return true;
	}
	return false; // No matches
}


/// @brief Has a Character?
bool
Cstring::has( char const c ) const
{
	for ( size_type i = 0; i < std::strlen( str_ ); ++i ) {
		if ( str_[ i ] == c ) return true;
	}
	return false; // No matches
}


/// @brief Length Space-Trimmed
Cstring::size_type
Cstring::len_trim() const
{
	for ( size_type i = std::strlen( str_ ); i > 0; --i ) {
		if ( str_[ i - 1 ] != SPACE ) return i;
	}
	return 0;
}


/// @brief Length Whitespace-Trimmed
Cstring::size_type
Cstring::len_trim_whitespace() const
{
	for ( size_type i = std::strlen( str_ ); i > 0; --i ) {
		char const c( str_[ i - 1 ] );
		if ( ( c != SPACE ) && ( c != TAB ) && ( c != NULLC ) ) return i;
	}
	return 0;
}


/// @brief Find First Occurrence of a Character
Cstring::size_type
Cstring::find( char const c ) const
{
	for ( size_type i = 0; i < std::strlen( str_ ); ++i ) {
		if ( str_[ i ] == c ) return i;
	}
	return npos; // No matches
}


/// @brief Find Last Occurrence of a Character
Cstring::size_type
Cstring::find_last( char const c ) const
{
	for ( size_type i = std::strlen( str_ ); i > 0; --i ) {
		if ( str_[ i - 1 ] == c ) return i;
	}
	return npos; // No matches
}


/// @brief Lowercase
Cstring &
Cstring::lowercase()
{
	for ( size_type i = 0; i < std::strlen( str_ ); ++i ) {
		str_[ i ] = std::tolower( str_[ i ] );
	}
	return *this;
}


/// @brief Uppercase
Cstring &
Cstring::uppercase()
{
	for ( size_type i = 0; i < std::strlen( str_ ); ++i ) {
		str_[ i ] = std::toupper( str_[ i ] );
	}
	return *this;
}


/// @brief Left Justify
Cstring &
Cstring::left_justify()
{
	size_type const len( std::strlen( str_ ) );
	for ( size_type i = 0; i < len; ++i ) {
		if ( str_[ i ] != SPACE ) {
			if ( i > 0 ) {
				std::memmove( str_, str_ + i, len - i );
				std::memset( str_ + len - i, SPACE, i );
			}
			return *this;
		}
	}
	return *this;
}


/// @brief Right Justify
Cstring &
Cstring::right_justify()
{
	size_type const len( std::strlen( str_ ) );
	for ( size_type i = len; i > 0; --i ) {
		if ( str_[ i - 1 ] != SPACE ) {
			if ( i < len ) {
				std::memmove( str_ + len - i, str_, i );
				std::memset( str_, SPACE, len - i );
			}
			return *this;
		}
	}
	return *this;
}


/// @brief Center
Cstring &
Cstring::center()
{
	left_justify();
	size_type const len_t( len_trim() );
	size_type const pad( ( length() - len_t ) / 2 );
	if ( pad > 0 ) {
		std::memmove( str_ + pad, str_, len_t );
		std::memset( str_, SPACE, pad );
	}
	return *this;
}


/// @brief Compress Out Whitespace
Cstring &
Cstring::compress()
{
	size_type const len( std::strlen( str_ ) );
	size_type j = 0;
	for ( size_type i = 0; i < len; ++i ) {
		char const c( str_[ i ] );
		if ( ( c != SPACE ) && ( c != TAB ) && ( c != NULLC ) ) str_[ j++ ] = c;
	}
	if ( j < len ) std::memset( str_ + j, SPACE, len - j );
	return *this;
}


/// @brief Cstring == Cstring Case-Insensitively?
bool
equali( Cstring const & s, Cstring const & t )
{
	using size_type = Cstring::size_type;
	size_type const s_len( s.length() );
	if ( s_len != t.length() ) {
		return false;
	} else {
		for ( size_type i = 0; i < s_len; ++i ) {
			if ( std::tolower( s.str_[ i ] ) != std::tolower( t.str_[ i ] ) ) return false;
		}
		return true;
	}
}


/// @brief Cstring == cstring Case-Insensitively?
bool
equali( Cstring const & s, c_cstring const t )
{
	using size_type = Cstring::size_type;
	size_type const s_len( s.length() );
	if ( s_len != std::strlen( t ) ) {
		return false;
	} else {
		for ( size_type i = 0; i < s_len; ++i ) {
			if ( std::tolower( s.str_[ i ] ) != std::tolower( t[ i ] ) ) return false;
		}
		return true;
	}
}


/// @brief cstring == Cstring Case-Insensitively?
bool
equali( c_cstring const s, Cstring const & t )
{
	using size_type = Cstring::size_type;
	size_type const s_len( std::strlen( s ) );
	if ( s_len != t.length() ) {
		return false;
	} else {
		for ( size_type i = 0; i < s_len; ++i ) {
			if ( std::tolower( s[ i ] ) != std::tolower( t.str_[ i ] ) ) return false;
		}
		return true;
	}
}


/// @brief Cstring == std::string Case-Insensitively?
bool
equali( Cstring const & s, std::string const & t )
{
	using size_type = Cstring::size_type;
	size_type const s_len( s.length() );
	if ( s_len != t.length() ) {
		return false;
	} else {
		for ( size_type i = 0; i < s_len; ++i ) {
			if ( std::tolower( s.str_[ i ] ) != std::tolower( t[ i ] ) ) return false;
		}
		return true;
	}
}


/// @brief std::string == Cstring Case-Insensitively?
bool
equali( std::string const & s, Cstring const & t )
{
	using size_type = Cstring::size_type;
	size_type const s_len( s.length() );
	if ( s_len != t.length() ) {
		return false;
	} else {
		for ( size_type i = 0; i < s_len; ++i ) {
			if ( std::tolower( s[ i ] ) != std::tolower( t.str_[ i ] ) ) return false;
		}
		return true;
	}
}


/// @brief Cstring == char Case-Insensitively?
bool
equali( Cstring const & s, char const c )
{
	return ( ( s.length() == 1 ) && ( std::tolower( s.str_[ 0 ] ) == std::tolower( c ) ) );
}


/// @brief char == Cstring Case-Insensitively?
bool
equali( char const c, Cstring const & s )
{
	return ( ( s.length() == 1 ) && ( std::tolower( s.str_[ 0 ] ) == std::tolower( c ) ) );
}


/// @brief Output to Stream
std::ostream &
operator <<( std::ostream & stream, Cstring const & s )
{
	for ( Cstring::size_type i = 0; i < std::strlen( s.str_ ); ++i ) {
		stream << s.str_[ i ];
	}
	return stream;
}


/// @brief Input from Stream
std::istream &
operator >>( std::istream & stream, Cstring & s )
{
	std::string buffer;
	stream >> buffer;
	s = buffer;
	return stream;
}


// Cstring


// Static Data Member Definitions

#ifndef _MSC_VER // Microsoft Visual C++ doesn't need or like these

Cstring::size_type const Cstring::npos;

#endif // _MSC_VER


} // namespace ObjexxFCL
