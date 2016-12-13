context("Code is high quality and lint free")

test_that("Code Lint", {
    skip_if_not_installed("lintr")
    my_linters <- list(
        absolute_paths_linter=lintr::absolute_paths_linter # Should not have absolute path in the code.
        # , assignment_linter=lintr::assignment_linter # Use <-, not =, for assignment.
        # , closed_curly_linter=lintr::closed_curly_linter # Closing curly-braces should always be on their own line, unless it's followed by an else.
        # , commas_linter=lintr::commas_linter # Commas should always have a space after.
        # , commented_code_linter=lintr::commented_code_linter 
        # , infix_spaces_linter=lintr::infix_spaces_linter # Put spaces around all infix operators.
        # , line_length_linter=lintr::line_length_linter(80) # line length of both comments and code is less than length.
        # , no_tab_linter=lintr::no_tab_linter
        # , object_usage_linter=lintr::object_usage_linter
        # , snake_case_linter=lintr::snake_case_linter
        # , multiple_dots_linter=lintr::multiple_dots_linter
        # , object_length_linter=lintr::object_length_linter # Objects do are not very long.not have.multiple.dots.
        # , open_curly_linter=lintr::open_curly_linter # Opening curly braces are never on their own line and are always followed by a newline.
        # , single_quotes_linter=lintr::single_quotes_linter
        # , spaces_inside_linter=lintr::spaces_inside_linter
        # , spaces_left_parentheses_linter=lintr::spaces_left_parentheses_linter
        # , trailing_blank_lines_linter=lintr::trailing_blank_lines_linter
        # , trailing_whitespace_linter=lintr::trailing_whitespace_linter
    )
    lintr::expect_lint_free(linters=my_linters) # uncomment this if you want to check code quality
})
