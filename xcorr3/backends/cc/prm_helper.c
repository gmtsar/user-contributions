#include "xcorr2.h"
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/* Keep source locations with values so conversion errors identify the field. */
struct prm_value {
    char *text;
    size_t line;
};

static void prm_value_free(void *data) {
    struct prm_value *value = data;
    if (value != NULL) {
        free(value->text);
        free(value);
    }
}

/* Parsing and conversion report failure; the caller owns process termination
 * and can therefore release every open PRM and other acquired resources. */
static bool prm_error(struct prm_handler *handler, size_t line_number,
                      const char *format, ...) {
    va_list args;
    if (handler == NULL)
        return false;
    snprintf(handler->error, sizeof(handler->error), "PRM %.240s",
             handler->filename != NULL ? handler->filename : "(unknown)");
    size_t used = strlen(handler->error);
    if (line_number != 0)
        snprintf(handler->error + used, sizeof(handler->error) - used,
                 ":%zu", line_number);
    used = strlen(handler->error);
    snprintf(handler->error + used, sizeof(handler->error) - used, ": ");
    used = strlen(handler->error);
    va_start(args, format);
    vsnprintf(handler->error + used, sizeof(handler->error) - used, format, args);
    va_end(args);
    return false;
}

/* Both scans stay within the string, including an empty/whitespace-only one. */
static char *trim(char *text) {
    while (isspace((unsigned char)*text))
        ++text;
    char *end = text + strlen(text);
    while (end > text && isspace((unsigned char)end[-1]))
        --end;
    *end = '\0';
    return text;
}

bool prm_open(struct prm_handler *handler, const char *fname) {
    if (handler == NULL)
        return false;
    memset(handler, 0, sizeof(*handler));
    if (fname == NULL || *fname == '\0')
        return prm_error(handler, 0, "missing filename");
    handler->filename = strdup(fname);
    if (handler->filename == NULL)
        return prm_error(handler, 0, "cannot allocate filename");

    FILE *fin = fopen(fname, "r");
    char *line = NULL;
    size_t capacity = 0, line_number = 0;
    ssize_t length;
    if (fin == NULL) {
        prm_error(handler, 0, "cannot open: %s", strerror(errno));
        goto fail;
    }
    handler->entry = g_hash_table_new_full(g_str_hash, g_str_equal,
                                          free, prm_value_free);
    while ((length = getline(&line, &capacity, fin)) >= 0) {
        ++line_number;
        if (memchr(line, '\0', (size_t)length) != NULL) {
            prm_error(handler, line_number, "embedded NUL byte");
            goto fail;
        }
        char *key = trim(line);
        if (*key == '\0' || *key == '#')
            continue;
        char *separator = strchr(key, '=');
        if (separator == NULL) {
            prm_error(handler, line_number, "expected key = value assignment");
            goto fail;
        }
        *separator = '\0';
        key = trim(key);
        if (*key == '\0') {
            prm_error(handler, line_number, "empty key");
            goto fail;
        }
        for (const char *p = key; *p != '\0'; ++p)
            if (isspace((unsigned char)*p)) {
                prm_error(handler, line_number, "whitespace inside key");
                goto fail;
            }
        char *text = trim(separator + 1);
        char *owned_key = strdup(key);
        struct prm_value *value = calloc(1, sizeof(*value));
        if (owned_key == NULL || value == NULL) {
            free(owned_key);
            free(value);
            prm_error(handler, line_number, "cannot allocate field");
            goto fail;
        }
        value->text = strdup(text);
        value->line = line_number;
        if (value->text == NULL) {
            free(owned_key);
            prm_value_free(value);
            prm_error(handler, line_number, "cannot allocate value");
            goto fail;
        }
        /* GMTSAR tools can append updated fields. Preserve the previous
         * last-assignment-wins behavior, including an explicitly empty value. */
        g_hash_table_replace(handler->entry, owned_key, value);
    }
    if (ferror(fin) || !feof(fin)) {
        prm_error(handler, line_number, "read failed: %s",
                  errno != 0 ? strerror(errno) : "I/O error");
        goto fail;
    }
    free(line);
    line = NULL;
    int close_status = fclose(fin);
    fin = NULL;
    if (close_status != 0) {
        prm_error(handler, line_number, "close failed: %s", strerror(errno));
        goto fail;
    }
    return true;

fail:
    free(line);
    if (fin != NULL)
        fclose(fin);
    prm_close(handler);
    return false;
}

void prm_close(struct prm_handler *handler) {
    if (handler == NULL)
        return;
    if (handler->entry != NULL)
        g_hash_table_destroy(handler->entry);
    free(handler->filename);
    handler->entry = NULL;
    handler->filename = NULL;
}

static const struct prm_value *required_value(struct prm_handler *handler,
                                              const char *key) {
    if (handler == NULL || handler->entry == NULL || key == NULL || *key == '\0') {
        prm_error(handler, 0, "invalid required-field lookup");
        return NULL;
    }
    const struct prm_value *value = g_hash_table_lookup(handler->entry, key);
    if (value == NULL) {
        prm_error(handler, 0, "missing required field '%s'", key);
        return NULL;
    }
    if (*value->text == '\0') {
        prm_error(handler, value->line, "empty required field '%s'", key);
        return NULL;
    }
    handler->error[0] = '\0';
    return value;
}

bool prm_get_str(struct prm_handler *handler, const char *key, const char **out) {
    if (out == NULL)
        return prm_error(handler, 0, "missing string output pointer");
    const struct prm_value *value = required_value(handler, key);
    if (value == NULL)
        return false;
    *out = value->text;
    return true;
}

bool prm_get_int(struct prm_handler *handler, const char *key, int *out) {
    if (out == NULL)
        return prm_error(handler, 0, "missing integer output pointer");
    const struct prm_value *value = required_value(handler, key);
    if (value == NULL)
        return false;
    char *end;
    errno = 0;
    long result = strtol(value->text, &end, 10);
    if (end == value->text || *end != '\0' || errno == ERANGE ||
        result < INT_MIN || result > INT_MAX)
        return prm_error(handler, value->line,
                         "field '%s' must be a decimal integer in [%d, %d]",
                         key, INT_MIN, INT_MAX);
    *out = (int)result;
    return true;
}

bool prm_get_f64(struct prm_handler *handler, const char *key, double *out) {
    if (out == NULL)
        return prm_error(handler, 0, "missing number output pointer");
    const struct prm_value *value = required_value(handler, key);
    if (value == NULL)
        return false;
    char *end;
    errno = 0;
    double result = strtod(value->text, &end);
    if (end == value->text || *end != '\0' || errno == ERANGE || !isfinite(result))
        return prm_error(handler, value->line,
                         "field '%s' must be a finite, representable number", key);
    *out = result;
    return true;
}
