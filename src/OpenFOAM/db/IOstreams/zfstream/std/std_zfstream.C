#include "std_zfstream.H"
#include "IOstreams.H"

// Construct a gzfilebuf object.
// Allocate memory for 'get' buffer and zero all buffer pointers.
gzfilebuf::gzfilebuf()
:
    std::streambuf(),
    file(NULL),
    mode(std::ios_base::openmode(0)),
    own_file_descriptor(false),
    ibuf_size(page_size/sizeof(char)),
    ibuffer(new char[ibuf_size])
{
    // Null get and set pointers.
    setg(0, 0, 0);
    setp(ibuffer, ibuffer + ibuf_size);
}


gzfilebuf::~gzfilebuf()
{
    sync();

    if (own_file_descriptor)
    {
        close();
    }

    delete[] ibuffer;
}


std::string gzfilebuf::modestr(std::ios_base::openmode m)
{
    std::string modeString;

    if (m & std::ios_base::in)
    {
        mode = std::ios_base::in;
        modeString += 'r';
    }
    else if (m & std::ios_base::app)
    {
        mode = std::ios_base::app;
        modeString += 'a';
    }
    else
    {
        mode = std::ios_base::out;
        modeString += 'w';
    }

    if (m & std::ios_base::binary)
    {
        mode |= std::ios_base::binary;
        modeString += 'b';
    }

    // Hard code the default compression level
    if (m & (std::ios_base::out | std::ios_base::app))
    {
        modeString += '1';
    }

    return modeString;
}


gzfilebuf* gzfilebuf::open
(
    const char *name,
    std::ios_base::openmode m
)
{
    if (is_open())
    {
        return NULL;
    }

    if ((file = gzopen(name, modestr(m).c_str())) == NULL)
    {
        return NULL;
    }

    own_file_descriptor = true;

    return this;
}


gzfilebuf* gzfilebuf::attach
(
    int file_descriptor,
    std::ios_base::openmode m
)
{
    if (is_open())
    {
        return NULL;
    }

    if ((file = gzdopen(file_descriptor, modestr(m).c_str())) == NULL)
    {
        return NULL;
    }

    own_file_descriptor = false;

    return this;
}


gzfilebuf* gzfilebuf::close()
{
    if (is_open())
    {
        sync();
        gzclose(file);
        file = NULL;
    }

    return this;
}


int gzfilebuf::setcompressionlevel(int comp_level)
{
    return gzsetparams(file, comp_level, -2);
}


int gzfilebuf::setcompressionstrategy(int comp_strategy)
{
    return gzsetparams(file, -2, comp_strategy);
}


gzfilebuf::pos_type gzfilebuf::seekoff
(
    gzfilebuf::off_type,
    std::ios_base::seekdir,
    std::ios_base::openmode
)
{
    return pos_type(traits_type::eof());
}


gzfilebuf::int_type gzfilebuf::overflow(int_type c)
{
    // Error if the file not open for writing
    if (!is_open() || !(mode & std::ios_base::out))
    {
        return traits_type::eof();
    }

    if (pptr() != 0)
    {
        if (pptr() > pbase())
        {
            if (flushbuf() == traits_type::eof())
            {
                return traits_type::eof();
            }
        }

        if (c != traits_type::eof())
        {
            *pptr() = c;
            pbump(1);
        }
    }

    return 0;
}


int gzfilebuf::sync()
{
    if (!is_open())
    {
        return traits_type::eof();
    }

    if (pptr() != 0 && pptr() > pbase())
    {
        return flushbuf();
    }

    return 0;
}


int gzfilebuf::flushbuf()
{
    char* q = pbase();
    int n = pptr() - q;

    if (gzwrite(file, q, n) < n)
    {
        // disable put area
        setp(0, 0);

        return traits_type::eof();
    }

    // Set the output (put) pointers
    setp(ibuffer, ibuffer + ibuf_size);

    return 0;
}


gzfilebuf::int_type gzfilebuf::underflow()
{
    // Error if the file not open for reading.
    if (!is_open() || !(mode & std::ios_base::in))
    {
        return traits_type::eof();
    }

    // If the input buffer is empty then try to fill it.
    if (gptr() != 0 && gptr() < egptr())
    {
        return (unsigned char) *gptr();
    }
    else
    {
        //Foam::Info<< int(gptr()) << " " << int(egptr()) << Foam::endl;

        return 
            fillbuf()
         == traits_type::eof() ? traits_type::eof() : (unsigned char) *gptr();
    }
}


// Load the input buffer from the underlying gz file.
// Returns number of characters read, or traits_type::eof().
//
int gzfilebuf::fillbuf()
{
    int t = gzread(file, ibuffer, ibuf_size);

    if (t <= 0)
    {
        // disable get area
        setg(0, 0, 0);
        return traits_type::eof();
    }

    // Set the input (get) pointers
    setg(ibuffer, ibuffer, ibuffer+t);

    return t;
}


gzifstream::gzifstream()
:
    std::istream(&buffer_)
{}


gzifstream::gzifstream(const char *name, std::ios_base::openmode m)
:
    std::istream(&buffer_)
{
    open(name, m);
}


gzifstream::gzifstream(int fd, std::ios_base::openmode m)
:
    std::istream(&buffer_)
{
    buffer_.attach(fd, m);
}


gzifstream::~gzifstream()
{}


void gzifstream::open(const char *name, std::ios_base::openmode m)
{
    if (!buffer_.open(name, m))
    {
        clear(std::ios_base::failbit | std::ios_base::badbit);
    }
    else
    {
        clear();
    }
}


void gzifstream::close()
{
    if (!buffer_.close())
    {
        clear(std::ios_base::failbit | std::ios_base::badbit);
    }
}


gzofstream::gzofstream()
:
    std::ostream(&buffer_)
{}


gzofstream::gzofstream(const char *name, std::ios_base::openmode m)
:
    std::ostream(&buffer_)
{
    open(name, m);
}


gzofstream::gzofstream(int fd, std::ios_base::openmode m)
:
    std::ostream(&buffer_)
{
    buffer_.attach(fd, m);
}


gzofstream::~gzofstream()
{}


void gzofstream::open(const char *name, std::ios_base::openmode m)
{
    if (!buffer_.open(name, m))
    {
        clear(std::ios_base::failbit | std::ios_base::badbit);
    }
    else
    {
        clear();
    }
}


void gzofstream::close()
{
    if (!buffer_.close())
    {
        clear(std::ios_base::failbit | std::ios_base::badbit);
    }
}
