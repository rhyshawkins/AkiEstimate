
#include "logging.hpp"

#include <time.h>

std::string logging::log::out("");
std::stringstream logging::log::tsbuffer;

logging::log::log()
{
}

bool
logging::log::set_output(const char *filename)
{
  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    return false;
  }

  fclose(fp);
  out = filename;

  return true;
}

std::string
logging::log::mkformatstring(const char *fmt, ...)
{
  static char *buffer = nullptr;
  static int buffer_size = -1;

  if (buffer == nullptr) {
    buffer_size = 512;
    buffer = new char[buffer_size];
  }

  va_list ap;
  int size;
  
  va_start(ap, fmt);
  size = vsnprintf(buffer, buffer_size, fmt, ap);
  while (size >= buffer_size) {
    delete [] buffer;
    buffer_size *= 2;
    buffer = new char[buffer_size];
    size = vsnprintf(buffer, buffer_size, fmt, ap);
  }
  va_end(ap);

  return std::string(buffer);
}

void
logging::log::vlog(const char *prefix,
		   const char *sourcefile,
		   const char *function,
		   int lineno,
		   const char *fmt,
		   va_list ap)
{
  FILE *fp;
  if (out.length() > 0) {
    fp = fopen(out.c_str(), "a");
    if (fp == NULL) {
      fprintf(stderr, "log::vlog: failed to open file %s\n", out.c_str());
      throw logging::fatalexception();
    }
  } else {
    fp = stderr;
  }

  fprintf(fp,
	  "%s:%s:%s:%s:%4d:",
	  timestamp(),
	  prefix,
	  sourcefile,
	  function,
	  lineno);
  vfprintf(fp, fmt, ap);
  fprintf(fp, "\n");

  if (out.length() > 0) {
    fclose(fp);
  }
}
		     
void
logging::log::vprintf(const char *fmt,
		      va_list ap)
{
  FILE *fp;
  if (out.length() > 0) {
    fp = fopen(out.c_str(), "a");
    if (fp == NULL) {
      fprintf(stderr, "log::vprintf: failed to open file %s\n", out.c_str());
      throw logging::fatalexception();
    }
  } else {
    fp = stderr;
  }

  vfprintf(fp, fmt, ap);

  if (out.length() > 0) {
    fclose(fp);
  }
}

void
logging::log::mark(const char *fmt, va_list ap)
{
  FILE *fp;
  if (out.length() > 0) {
    fp = fopen(out.c_str(), "a");
    if (fp == NULL) {
      fprintf(stderr, "log::mark: failed to open file %s\n", out.c_str());
      throw logging::fatalexception();
    }
  } else {
    fp = stderr;
  }

  fprintf(fp, "%s:", timestamp());
  vfprintf(fp, fmt, ap);
  fprintf(fp, "\n");

  if (out.length() > 0) {
    fclose(fp);
  }
}

const char *
logging::log::timestamp()
{
  const char *TIME_FORMAT = "%Y-%m-%d %H:%M:%S";
  static char buffer[256];
  time_t tmp;
  struct tm *t;

  tmp = time(NULL);
  t = localtime(&tmp);
  
  if (strftime(buffer, sizeof(buffer), TIME_FORMAT, t) == 0) {
    fprintf(stderr, "log::timestamp: failed to format time\n");
    return NULL;
  }

  return buffer;
}
