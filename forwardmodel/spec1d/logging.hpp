#pragma once
#ifndef logging_hpp
#define logging_hpp

#include <stdarg.h>
#include <stdio.h>

#include <iostream>
#include <sstream>

#define FATAL(fmt, ...) logging::fatal::log(__FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define ERROR(fmt, ...) logging::error::log(__FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define WARNING(fmt, ...) logging::warning::log(__FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define INFO(fmt, ...) logging::info::log(__FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)
#define DEBUG(fmt, ...) logging::debug::log(__FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)

#define STATUS(fmt, ...) logging::status::log(fmt, ##__VA_ARGS__)
#define MARK(fmt, ...) logging::status::mark(fmt, ##__VA_ARGS__)

namespace logging {
  
  class fatalexception : public std::exception {
  };
    
  class log {
  public:
    
    static bool set_output(const char *filename);

    static std::string mkformatstring(const char *fmt, ...);
    
  protected:
    
    static void vlog(const char *prefix,
		     const char *sourcefile,
		     const char *function,
		     int lineno,
		     const char *fmt,
		     va_list ap);

    static void vprintf(const char *fmt,
			va_list ap);

    static void mark(const char *fmt,
		     va_list ap);
		     
  private:
    log();

    static const char *timestamp();
    
    static std::string out;
    static std::stringstream tsbuffer;
  };

  class fatal : public log {
  public:
    
    static void log(const char *sourcefile,
		    const char *function,
		    int lineno,
		    const char *fmt,
		    ...)
    {
      va_list ap;

      va_start(ap, fmt);
      log::vlog("fatal", sourcefile, function, lineno, fmt, ap);
      va_end(ap);

      throw fatalexception();
    }
  };
  
  class error : public log {
  public:

    static void log(const char *sourcefile,
		    const char *function,
		    int lineno,
		    const char *fmt,
		    ...)
    {
      va_list ap;

      va_start(ap, fmt);
      log::vlog("error", sourcefile, function, lineno, fmt, ap);
      va_end(ap);
    }
  };

  class warning : public log {
  public:

    static void log(const char *sourcefile,
		    const char *function,
		    int lineno,
		    const char *fmt,
		    ...)
    {
      va_list ap;

      va_start(ap, fmt);
      log::vlog("warning", sourcefile, function, lineno, fmt, ap);
      va_end(ap);
    }
  };

  class info : public log {
  public:

    static void log(const char *sourcefile,
		    const char *function,
		    int lineno,
		    const char *fmt,
		    ...)
    {
      va_list ap;

      va_start(ap, fmt);
      log::vlog("info", sourcefile, function, lineno, fmt, ap);
      va_end(ap);
    }
  };

  class debug : public log {
  public:

    static void log(const char *sourcefile,
		    const char *function,
		    int lineno,
		    const char *fmt,
		    ...)
    {
      va_list ap;

      va_start(ap, fmt);
      log::vlog("debug", sourcefile, function, lineno, fmt, ap);
      va_end(ap);
    }
  };

  class status : public log {
  public:

    static void log(const char *fmt,
		    ...)
    {
      va_list ap;

      va_start(ap, fmt);
      log::vprintf(fmt, ap);
      va_end(ap);
    }

    static void mark(const char *fmt,
		    ...)
    {
      va_list ap;

      va_start(ap, fmt);
      log::mark(fmt, ap);
      va_end(ap);
    }
      
  };
  
  

};

#endif // logging
