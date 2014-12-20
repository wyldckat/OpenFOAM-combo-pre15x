#ifndef _zfstream_h
#define _zfstream_h

#include <iostream>
#include "zlib.h"

class gzfilebuf
:
    public std::streambuf
{
    gzFile file;
    std::ios_base::openmode mode;
    bool own_file_descriptor;

    // Get (input) buffer.
    int ibuf_size;
    char* ibuffer;

    enum
    {
        page_size = 4096
    };

    int flushbuf();
    int fillbuf();

    // Convert io_mode to "rwab" string.
    std::string modestr(std::ios_base::openmode io_mode);


    // Not defined
    gzfilebuf(const gzfilebuf&);
    void operator=(const gzfilebuf&);


protected:

    virtual int_type underflow();
    virtual int_type overflow(int_type c = traits_type::eof());


public:

    gzfilebuf();
    virtual ~gzfilebuf();

    gzfilebuf* open(const char* name, std::ios_base::openmode io_mode);
    gzfilebuf* attach(int file_descriptor, std::ios_base::openmode io_mode);
    gzfilebuf* close();

    int setcompressionlevel
    (
        int comp_level = Z_BEST_SPEED
        /*
        Z_NO_COMPRESSION = 0
        Z_BEST_SPEED = 1
        Z_BEST_COMPRESSION = 9
        Z_DEFAULT_COMPRESSION = 6
        */
    );

    int setcompressionstrategy
    (
        int comp_strategy = Z_DEFAULT_STRATEGY
        /*
        Z_FILTERED            1
        Z_HUFFMAN_ONLY        2
        Z_DEFAULT_STRATEGY    0
        */
    );

    bool is_open() const
    {
        return (file != NULL);
    }

    virtual pos_type seekoff
    (
        off_type off,
        std::ios_base::seekdir way,
        std::ios_base::openmode m = std::ios_base::in|std::ios_base::out
    );

    virtual int sync();
};


class gzifstream
:
    public std::istream
{
    gzfilebuf buffer_;

public:

    gzifstream();
    gzifstream(const char *name, std::ios_base::openmode io_mode = ios_base::in);
    gzifstream(int fd, std::ios_base::openmode io_mode = ios_base::in);

    virtual ~gzifstream();

    void open(const char *name, std::ios_base::openmode io_mode);
    void close();
};


class gzofstream
:
    public std::ostream
{

public:

    gzfilebuf buffer_;

    gzofstream();
    gzofstream(const char *name, std::ios_base::openmode io_mode = ios_base::out);
    gzofstream(int fd, std::ios_base::openmode io_mode = ios_base::out);

    virtual ~gzofstream();

    void open(const char *name, std::ios_base::openmode io_mode);
    void close();
};


template<class T>
class gzomanip;

template<class T>
gzofstream &operator<<(gzofstream &, const gzomanip<T> &);

template<class T>
class gzomanip
{
    friend gzofstream &operator<< <T>(gzofstream &, const gzomanip<T> &);

public:

    gzomanip(gzofstream &(*f)(gzofstream &, T), T v)
    :
        func(f),
        val(v)
    {}

private:

    gzofstream &(*func)(gzofstream &, T);
    T val;
};


template<class T> gzofstream &operator<<
(
    gzofstream &s,
    const gzomanip<T> &m
)
{
    return (*m.func)(s, m.val);
}

inline gzofstream& setcompressionlevel(gzofstream& s, int l)
{
    s.buffer_.setcompressionlevel(l);
    return s;
}

inline gzofstream& setcompressionstrategy(gzofstream& s, int l)
{
    s.buffer_.setcompressionstrategy(l);
    return s;
}

inline gzomanip<int> setcompressionlevel(int l)
{
    return gzomanip<int>(&setcompressionlevel, l);
}

inline gzomanip<int> setcompressionstrategy(int l)
{
    return gzomanip<int>(&setcompressionstrategy, l);
}

#endif
