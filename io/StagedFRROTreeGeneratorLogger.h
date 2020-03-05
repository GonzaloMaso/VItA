#ifndef STAGEDFRROTREEGENERATORLOGGER_H
#define STAGEDFRROTREEGENERATORLOGGER_H

class StagedFRROTreeGeneratorLogger 
{
    // Attibutes
    StagedFRROTreeGenerator *tree_generator_;
    FILE *file_;

    public:
        /**
         * Constructor.
         * @param fp Pointer to file to write to, in "w" mode.
         * @param tree_gen Tree generator.
         */
        StagedFRROTreeGeneratorLogger(FILE *fp, StagedFRROTreeGenerator* tree_gen);
        ~StagedFRROTreeGeneratorLogger();
        /**
         * Writes relevant parameters of @p tree_generator_ to @p file_.
         */
        void write();
};

#endif // STAGEDFRROTREEGENERATORLOGGER_H