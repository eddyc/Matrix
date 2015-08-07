//
//  Matrix.h
//  ModVoc
//
//  Created by Edward Costello on 23/08/2014.
//  Copyright (c) 2014 Edward Costello. All rights reserved.
//


#import <MacTypes.h>
#import <stdio.h>
#import <stdlib.h>
#import <hdf5.h>
#import <Block.h>
#import <complex.h>
#import "CommonDSP.h"

#ifdef __cplusplus
extern "C"
{
#endif
    
#pragma mark - Structure -
    
    typedef struct Matrix
    {
        size_t elementCount;
        const size_t allocatedElementCount;
        size_t rowCount;
        size_t columnCount;
        size_t rowStride;
        size_t columnStride;
        Float64 *data;
        Float64 *realData;
        Float64 *imagData;
        Float64 *tempData;
        Float64 *tempImagData;
        bool *isComplex;
        bool isSubmatrix;
        bool dataOnStack;
        char *label;
    } Matrix;
   
#pragma mark - Allocation / Deallocation -
    
    OVERLOADED Matrix *Matrix_new(size_t rowCount , size_t columnCount);
    OVERLOADED Matrix *Matrix_new(size_t rowCount , size_t columnCount, char label[]);
    OVERLOADED Matrix *Matrix_new(Matrix *input);
    OVERLOADED Matrix Matrix_new(Float64 data[], size_t rowCount, size_t columnCount);
    OVERLOADED Matrix Matrix_new(Float64 data[], size_t columnCount);
    void Matrix_delete(Matrix *self);
    void Matrix_scopedDelete(Matrix **self);
#define _Matrix __attribute__((unused)) __attribute__((cleanup(Matrix_scopedDelete))) Matrix
    
#pragma mark - Setters / Getters -
   
    typedef struct { size_t value; } RowCount;
    typedef struct { size_t value; } ColumnCount;
    typedef struct { size_t value; } RowOffset;
    typedef struct { size_t value; } ColumnOffset;
    typedef struct { size_t value; } SourceRowOffset;
    typedef struct { size_t value; } SourceColumnOffset;
    typedef struct { size_t value; } DestinationRowOffset;
    typedef struct { size_t value; } DestinationColumnOffset;
    
    OVERLOADED Float64 Matrix_getElement(Matrix *input, size_t element);
    OVERLOADED Float64 Matrix_getElement(Matrix *input, RowOffset rowOffset, ColumnOffset columnOffset);
    
    OVERLOADED void Matrix_setElement(Matrix *input, size_t index, Float64 value);
    OVERLOADED void Matrix_setElement(Matrix *input, RowOffset rowOffset, ColumnOffset columnOffset, Float64 value);
    
    OVERLOADED void Matrix_setComplex(Matrix *input);
    OVERLOADED void Matrix_setComplex(Matrix *input, bool state);

#pragma mark - Safety -
    
#ifdef DEBUG
#define _Matrix_checkDimensionEquality(count, ...) Matrix_checkDimensionEquality(count, __VA_ARGS__)
#define _Matrix_checkElementCount(self, elementCount) Matrix_checkElementCount(self, elementCount)
#define _Matrix_checkForVector(self) Matrix_checkForVector(self)
#define _Matrix_checkSubmatrixDimensions(input, rowCount, rowOffset, columnCount, columnOffset) Matrix_checkSubmatrixDimensions(input, rowCount, rowOffset, columnCount, columnOffset)
#define Matrix_getColumnCount(matrix) Matrix_getColumnCount(matrix)
#define _Matrix_checkForValidIndex(input, index) Matrix_checkForValidIndex(input, index)
#define _Matrix_checkIsEmpty(matrix) Matrix_checkIsEmpty(matrix)
#define _Matrix_checkCharacterEquality(a, b) Matrix_checkCharacterEquality(a, b)
#define _Matrix_checkIsComplex(matrix) Matrix_checkIsComplex(matrix)
#define _Matrix_checkIsNotComplex(matrix) Matrix_checkIsNotComplex(matrix)
#else
#define _Matrix_checkDimensionEquality(count, ...)
#define _Matrix_checkElementCount(self, elementCount)
#define _Matrix_checkForVector(self)
#define _Matrix_checkSubmatrixDimensions(input, rowCount, rowOffset, columnCount, columnOffset)
#define Matrix_getColumnCount(matrix) (matrix)->columnCount
#define _Matrix_checkForValidIndex(input, index)
#define _Matrix_checkIsEmpty(matrix)
#define _Matrix_checkCharacterEquality(a, b)
#define _Matrix_checkIsComplex(matrix)
#define _Matrix_checkIsNotComplex(matrix)
#endif
    
#define Matrix_getRow(matrix, row) (&(matrix)->data[(matrix)->rowStride * row])
#define Matrix_getImagRow(matrix, row) (&(matrix)->imagData[(matrix)->rowStride * row])
#define Matrix_isComplex(matrix) (matrix)->isComplex[0]
    
    typedef __block bool (^CompareElement)(Float64);
    typedef __block bool (^CompareIndex)(size_t);
    typedef __block Float64 (^ChangeElement)(Float64);
    
    OVERLOADED void Matrix_checkDimensionEquality(Matrix *inputA, Matrix *inputB);
    OVERLOADED void Matrix_checkDimensionEquality(size_t inputCount, ...);
    
    void Matrix_checkElementCount(Matrix *self, size_t elementCount);
    void Matrix_checkForVector(Matrix *self);
    void Matrix_checkSubmatrixDimensions(Matrix *input,
                                         size_t rowCount, size_t rowOffset,
                                         size_t columnCount, size_t columnOffset);
    Boolean Matrix_checkIsSubmatrix(Matrix *input);
    void Matrix_checkForValidIndex(Matrix *input, size_t index);
    void Matrix_checkIsEmpty(Matrix *self);
    
    void Matrix_checkIsComplex(Matrix *self);
#pragma mark - Data copy -
    
    OVERLOADED void Matrix_copy(Matrix source, Matrix destination);
    OVERLOADED void Matrix_copy(Matrix *source, Matrix *destination);
    OVERLOADED void Matrix_copyRawData(Matrix *self,
                                       size_t elementCount,
                                       Float64 *data);
    
    OVERLOADED void Matrix_copyRawComplexData(Matrix *self,
                                              size_t elementCount,
                                              Float64 *realData,
                                              Float64 *imagData);
    
    OVERLOADED void Matrix_subMatrixCopy(Matrix *source,
                                         ColumnCount columnCount,
                                         Matrix *destination);
    
    OVERLOADED void Matrix_subMatrixCopy(Matrix *source,
                                         RowCount rowCount,
                                         ColumnCount columnCount,
                                         Matrix *destination);
    
    OVERLOADED void Matrix_subMatrixCopy(Matrix *source,
                                         SourceRowOffset sourceStartRow,
                                         SourceColumnOffset sourceStartColumn,
                                         RowCount rowCount,
                                         ColumnCount columnCount,
                                         Matrix *destination,
                                         DestinationRowOffset destinationStartRow,
                                         DestinationColumnOffset destinationStartColumn);
    
#pragma mark - Submatrix View -
    
    Matrix Matrix_rowView(Matrix *input,
                          size_t rowOffset);
    Matrix Matrix_columnView(Matrix *input,
                          size_t columnOffset);
    OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                           RowCount rowCount);
    
    OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                           RowCount rowCount,
                                           RowOffset rowOffset);
    
    OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                           RowCount rowCount,
                                           ColumnCount columnCount);
    
    OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                           ColumnCount columnCount);
    OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                           ColumnCount columnCount,
                                           ColumnOffset columnOffset);
    OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                           ColumnOffset columnOffset);
    OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                           RowCount rowCount,
                                           ColumnCount columnCount,
                                           ColumnOffset columnOffset);
    
    OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                           RowCount rowCount,
                                           RowOffset rowOffset,
                                           ColumnCount columnCount);
    
    OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                           RowCount rowCount,
                                           RowOffset rowOffset,
                                           ColumnCount columnCount,
                                           ColumnOffset columnOffset);
    
    Matrix Matrix_tempView(Matrix *input);
    Matrix Matrix_realView(Matrix *input);
    Matrix Matrix_imagView(Matrix *input);
    Matrix Matrix_complexFlipView(Matrix *input);
    void Matrix_complexFlip(Matrix **input);

#pragma mark - Serialisation -
    
#define MatrixHDF(matrix) Matrix_saveHDF5File((matrix), "Desktop", "/"#matrix)
    
    OVERLOADED void Matrix_saveHDF5File(Matrix self, const char *path, const char *dataSet);
    OVERLOADED void Matrix_saveHDF5File(Matrix *self, const char *path, const char *dataSet);
    OVERLOADED Matrix *Matrix_openHDF5File(const char *path, const char *dataSet);
    void Matrix_addHDF5Dataset(Matrix *self, const char *path, const char *dataSet);
    
#pragma mark - Printing -
    
    OVERLOADED void Matrix_print(Matrix *input);
    OVERLOADED void Matrix_print(Matrix *input, const char *label);
    OVERLOADED void Matrix_print(Matrix input);
    OVERLOADED void Matrix_print(Matrix *input, CompareIndex comparison, size_t triggerIndex);
    
    void Matrix_printRange(const Matrix *const self,
                           size_t startRow,
                           size_t rowCount,
                           size_t startColumn,
                           size_t columnCount);
    
#pragma mark - Windowing functions -
    
    OVERLOADED void Matrix_hanningWindow(Matrix *input);
    OVERLOADED void Matrix_sineWindow(Matrix *input);
    OVERLOADED void Matrix_sineWindow(Matrix *input, Float64 startValue, Float64 increment);
    
#pragma mark - Generate -
    
    OVERLOADED void Matrix_reshape(Matrix *input, RowCount rows, ColumnCount columns);
    OVERLOADED void Matrix_reshape(Matrix *input, ColumnCount columns);
    OVERLOADED void Matrix_reshape(Matrix *input, RowCount rows);
    OVERLOADED void Matrix_reshape(Matrix *input, Matrix *shape);
    OVERLOADED void Matrix_transpose(Matrix *input);
    OVERLOADED void Matrix_transpose(Matrix *input, Matrix *output);
    
    OVERLOADED void Matrix_fill(Matrix input, Float64 scalar);
    OVERLOADED void Matrix_fill(Matrix *input, Float64 scalar);
    
    void Matrix_clear(const Matrix *const input);
    void Matrix_clearReal(const Matrix *const input);
    
    OVERLOADED void Matrix_ramp(Matrix input, Float64 startValue, Float64 increment);
    OVERLOADED void Matrix_ramp(Matrix *input, Float64 startValue, Float64 increment);
    OVERLOADED Matrix *Matrix_ramp(Float64 startValue, Float64 endValue, Float64 increment);
    
    OVERLOADED void Matrix_fillElements(Matrix *input, CompareElement comparison, Matrix *result);
    OVERLOADED void Matrix_fillElements(Matrix *input, CompareIndex comparison, Matrix *result);
    OVERLOADED void Matrix_fillElements(Matrix *input, CompareElement comparison, Matrix *comparisonData, Matrix *result);
    
    void Matrix_fillIndices(Matrix *input, CompareElement comparison, Matrix *result);
    OVERLOADED void Matrix_setElements(Matrix *input, CompareElement comparison, Float64 value, Matrix *result);
    OVERLOADED void Matrix_setElements(Matrix *input, CompareElement comparison,  Float64 (^value)(Float64 element), Matrix *result);
    void Matrix_insert(Matrix *input, Float64 value, size_t index);
    
    OVERLOADED void Matrix_unique(Matrix *input);
    OVERLOADED void Matrix_unique(Matrix *input, Matrix *output);
    
    OVERLOADED void Matrix_complexConjugate(Matrix *input, Matrix *output);
    
    void Matrix_gaussian(Matrix *input, Float64 mean, Float64 variance);
    void Matrix_diagonal(Matrix *input, Matrix *diagonal);
    void Matrix_identity(Matrix *input);
#pragma mark - Arithmetic functions -
    
    OVERLOADED void Matrix_add(Matrix input, Float64 scalar);
    OVERLOADED void Matrix_add(Matrix *input, Float64 scalar);
    OVERLOADED void Matrix_add(Matrix *input, Float64 scalar, Matrix *output);
    OVERLOADED void Matrix_add(Matrix *inputA,
                               Matrix *inputB,
                               Matrix *result);
    
    OVERLOADED void Matrix_subtract(Matrix *input,
                                    Matrix *negated,
                                    Matrix *result);
    OVERLOADED void Matrix_negate(Matrix *input);
    OVERLOADED void Matrix_negate(Matrix *input, Matrix *output);
    OVERLOADED void Matrix_multiply(Matrix input, Float64 scalar);
    OVERLOADED void Matrix_multiply(Matrix *input, Float64 scalar);
    OVERLOADED void Matrix_multiply(Matrix *input, Float64 scalar, Matrix *output);
    OVERLOADED void Matrix_multiply(Matrix *input,
                                    double complex scalar);
    OVERLOADED void Matrix_multiply(Matrix *input,
                                    double complex scalar,
                                    Matrix *output);
    OVERLOADED void Matrix_multiply(Matrix *inputA,
                                    Matrix *inputB,
                                    Matrix *result);
    void Matrix_dotProduct(Matrix *inputA, Matrix *inputB, Matrix *output);
    OVERLOADED void Matrix_divide(Matrix input, Float64 scalar);
    OVERLOADED void Matrix_divide(Matrix *input, Float64 scalar);
    OVERLOADED void Matrix_divide(Matrix *input, Float64 scalar, Matrix *output);
    OVERLOADED void Matrix_divide(Matrix *inputA,
                                  Matrix *inputB,
                                  Matrix *result);
    
    OVERLOADED Float64 Matrix_sum(Matrix *input, size_t startIndex, size_t size);
    OVERLOADED Float64 Matrix_sum(Matrix *input);
    OVERLOADED Float64 Matrix_sum(Matrix *input, Matrix *indices);
    void Matrix_sumRows(Matrix *input, Matrix *output);
    OVERLOADED void Matrix_difference(Matrix *input, Matrix *output);
    OVERLOADED void Matrix_difference(Matrix *inputA, Matrix *inputB, Matrix *output);
    Matrix *Matrix_newDifference(Matrix *inputA, Matrix *inputB);

    OVERLOADED void Matrix_floor(Matrix *input);
    OVERLOADED void Matrix_floor(Matrix *input, Matrix *output);
    OVERLOADED void Matrix_ceiling(Matrix *input);
    OVERLOADED void Matrix_ceiling(Matrix *input, Matrix *output);
    OVERLOADED void Matrix_cummulativeSum(Matrix *input);
    OVERLOADED void Matrix_cummulativeSum(Matrix *input, Matrix *output);
#pragma mark - Algebraic functions -
    
    OVERLOADED void Matrix_square(Matrix *input);
    OVERLOADED void Matrix_square(Matrix *input, Matrix *output);
    Float64 Matrix_sumOfSquares(Matrix *input);
    
    OVERLOADED void Matrix_squareRoot(Matrix *input);
    OVERLOADED void Matrix_squareRoot(Matrix *input, Matrix *output);
    void Matrix_modulus(Matrix *input, Matrix *output);
    
    OVERLOADED void Matrix_power(Matrix *input, Matrix *power, Matrix *output);
    OVERLOADED void Matrix_log10(Matrix *input);
    OVERLOADED void Matrix_log10(Matrix *input, Matrix *output);
    
#pragma mark - Trigonometric functions -
    
    void Matrix_angle(Matrix *input, Matrix *result);
    void Matrix_sine(Matrix *input, Matrix *result);
    void Matrix_cosine(Matrix *input, Matrix *result);
    void Matrix_complexExponent(Matrix *input, size_t frameNumber);
    void Matrix_complexExponentBlankReal(Matrix *input);
    void Matrix_magnitudeSquared(Matrix *input, Matrix *output);
    
#pragma mark - Information -
    
    void Matrix_minimum(Matrix *input, Float64 *minimum, size_t *index);
    void Matrix_maximum(Matrix *input, Float64 *maximum, size_t *index);
    void Matrix_domain(Matrix *input, Float64 domain[2]);
    size_t Matrix_countElements(Matrix *input, bool (^comparison)(Float64));
    
#pragma mark - Other functions -
    
    OVERLOADED void Matrix_interpolate(Matrix *input, Matrix *indices, Matrix *output);
    OVERLOADED void Matrix_interpolate(Matrix *currentIndices, Matrix *currentValues,
                                       Matrix *sortIndices, Matrix *sortedCurrentIndices, Matrix *sortedCurrentValues,
                                       Matrix *newIndices, Matrix *output);
    OVERLOADED void Matrix_fastInterpolate(Matrix *currentIndices, Matrix *currentValues,
                                           Matrix *newIndices, Matrix *output);
    OVERLOADED void Matrix_getSortIndices(Matrix *input, Matrix *indices, bool forward);
    OVERLOADED void Matrix_sortUsingIndices(Matrix *input, Matrix *indices);
    OVERLOADED void Matrix_sortUsingIndices(Matrix *input, Matrix *indices, Matrix *output);
    OVERLOADED void Matrix_sort(Matrix *input, Boolean forward);
    OVERLOADED void Matrix_sort(Matrix *input, Matrix *output, Boolean forward);
   
    
    OVERLOADED void Matrix_detrend(Matrix *input);
    OVERLOADED void Matrix_detrend(Matrix *input, Matrix *output);
    OVERLOADED void Matrix_round(Matrix *input, Matrix *output);
    OVERLOADED void Matrix_round(Matrix *input);
    
    Float64 Matrix_sumIndexedElements(Matrix *input, Matrix *indices);
    
    OVERLOADED void Matrix_absolute(Matrix *input);
    OVERLOADED void Matrix_absolute(Matrix *input, Matrix *output);
    
    OVERLOADED void Matrix_reverse(Matrix *input);
    OVERLOADED void Matrix_reverse(Matrix *input, Matrix *output);
    
#ifdef __cplusplus
}
#endif
