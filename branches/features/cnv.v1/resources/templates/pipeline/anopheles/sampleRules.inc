[%- sample = manifest.get_sample(sample_id) -%]

$(ALIGNMENTS_PATH)/[% sample_id %].bam : [% INCLUDE requirements %]
[% IF sample.dataSources.size >= 1 -%]
	$(PGV) merge $(addprefix -i , $^) -o $@
[%- ELSE -%]
	@echo missing alignment $@
	@exit 1
[%- END -%]

[%- IF sample.has_inout_data_source -%]
.PHONY: $(ALIGNMENTS_PATH)/[% sample_id %].bam
[%- END -%]

[%- FOREACH source IN sample.dataSources %]

$(TMP_DIR)/[% source.file_basename %]: $(TMP_DIR)/[% source.file_basename %].pre
	$(PGV) align -i $< -o $@

.SECONDARY: $(TMP_DIR)/[% source.file_basename %].pre

[%- END ~%]


[%~ BLOCK requirements -%]
[%- IF sample.has_inout_data_source -%]
[%- FOREACH data_source IN sample.dataSources -%]
| $(TMP_DIR)/[% data_source.file_basename %][% ' ' %]
[%- END -%]
[%- END -%]
[%- END -%]

