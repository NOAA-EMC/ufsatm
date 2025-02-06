## Description

(Instructions: this, and all subsequent sections of text should be removed and filled in as appropriate.)
Provide a detailed description of what this PR does.  
What bug does it fix, or what feature does it add?  
Is a change of answers expected from this PR?  



### Issue(s) addressed

Link the issues to be closed with this PR, whether in this repository, or in another repository.
(Remember, issues should always be created before starting work on a PR branch!)
- fixes #<issue_number>
- fixes noaa-emc/fv3atm/issues/<issue_number>



## Testing

How were these changes tested?  
What compilers / HPCs was it tested with?  
Are the changes covered by regression tests? (If not, why? Do new tests need to be added?)  
Have the ufs-weather-model regression test been run? On what platform?  
- Will the code updates change regression test baseline? If yes, why? Please show the baseline directory below.
- Please commit the regression test log files in your ufs-weather-model branch


## Dependencies

If testing this branch requires non-default branches in other repositories, list them.
Those branches should have matching names (ideally)

Do PRs in upstream repositories need to be merged first?
If so add the "waiting for other repos" label and list the upstream PRs
- waiting on noaa-emc/nems/pull/<pr_number>
- waiting on noaa-emc/fv3atm/pull/<pr_number>

# Requirements before merging
- [ ] All new code in this PR is tested by at least one unit test
- [ ] All new code in this PR includes Doxygen documentation
- [ ] All new code in this PR does not add new compilation warnings (check CI output)
