#ifndef STAGEDFRROTREEGENERATORLOGGER_H
#define STAGEDFRROTREEGENERATORLOGGER_H

class StagedFRROTreeGeneratorLogger 
{
    // Attibutes
    StagedFRROTreeGenerator *treeGenerator;
    FILE *file;

    public:
        /**
         * Constructor.
         * @param fp Pointer to file to write to, in "w" mode.
         * @param treeGen Tree generator.
         */
        StagedFRROTreeGeneratorLogger(FILE *fp, StagedFRROTreeGenerator* treeGen);
        ~StagedFRROTreeGeneratorLogger();
        /**
         * Writes relevant parameters of @p treeGenerator_ to @p file.
         */
        void write();
};

#endif // STAGEDFRROTREEGENERATORLOGGER_H