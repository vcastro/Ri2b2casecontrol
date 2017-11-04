/****** Object:  StoredProcedure [dbo].[usp_build_analysistables]    Script Date: 11/3/2017 3:18:12 PM ******/
DROP PROCEDURE [dbo].[usp_build_analysistables]

--stored procedure to write out analysis files for a patient_set to run case-control analysis
CREATE procedure [dbo].[usp_build_analysistables] 
(
	@outcome_concept_paths VARCHAR(MAX),
	@outcome_name VARCHAR(1000),
	@exposure_concept_paths VARCHAR(MAX),
	@visit_concept_paths VARCHAR(MAX),
	@patientset_instance_id INT = NULL
)
as

declare @parsed_outcome_concept_paths VARCHAR(MAX)
declare @parsed_exposure_concept_paths VARCHAR(MAX)
declare @parsed_visit_concept_paths VARCHAR(MAX)

--set @c_fullname_list = '\i2b2metadata\Diagnosis_ICD10\(M00-M99) Dise~6mvn\(M80-M85) Diso~xkfm\(M81) Osteopor~zztr\|\i2b2metadata\Medications_RxNorm\MRX\(N0000029177)~amx9\(N0000029178)~idib\%|\i2b2metadata\Medications_RxNorm\MRX\(N0000029177)~amx9\(N0000029178)~idib\(N0000029312)~n6af\(5492)~x6zg\|\i2b2metadata\Medications_RxNorm\MRX\(N0000029177)~amx9\(N0000029178)~idib\(N0000029312)~n6af\(8640)~0kua\'
--set @patientset_instance_id = 61

set @parsed_outcome_concept_paths	= '(''' + replace(@outcome_concept_paths, '|', ''', ''') + ''')'
set @parsed_exposure_concept_paths	= '(''' + replace(@exposure_concept_paths, '|', ''', ''') + ''')'
set @parsed_visit_concept_paths		= '(''' + replace(@visit_concept_paths, '|', ''', ''') + ''')'

declare @exec_string NVARCHAR(MAX)
declare @Rcd INT


-- output demographics
IF @patientset_instance_id IS NULL
BEGIN
	select p.patient_num, cast(p.birth_date as DATE) birth_date, cast(p.death_date as DATE) death_date, p.RACE_CD, p.SEX_CD, p.ETHNICITY_CD
	from patient_dimension p
END
ELSE
BEGIN
	select p.patient_num, cast(p.birth_date as DATE) birth_date, cast(p.death_date as DATE) death_date, p.RACE_CD, p.SEX_CD, p.ETHNICITY_CD
	from patient_dimension p, QT_PATIENT_SET_COLLECTION ps
	where p.PATIENT_NUM = ps.PATIENT_NUM and ps.RESULT_INSTANCE_ID = @patientset_instance_id
END

-- output dictionary
select @exec_string = 
	'select concept_cd, name_char
		from concept_dimension
		where concept_path in '+@parsed_exposure_concept_paths +' 
	union
	select ''V'', ''Visits''
	union
	select ''O'', '''+@outcome_name+''''

	EXEC @Rcd = sp_executesql @exec_string	  

	IF @Rcd <> 0 OR @Rcd is NULL
	begin
		print 'Failed in building dictionary'
		return isnull (@Rcd, 99) -- return error code to the calling program
	end	

-- output analysis table
select @exec_string =
	'select patient_num, isnull(c.concept_cd, left(c.name_char,50)) concept_cd, cast(start_date as DATE) concept_date, row_number() over (partition by patient_num, isnull(c.concept_cd, left(c.name_char,50)) order by cast(start_date as DATE)) concept_order
	from observation_fact f, concept_dimension c, concept_dimension c2
	where f.CONCEPT_CD = c2.CONCEPT_CD and c2.CONCEPT_PATH like c.CONCEPT_PATH + ''%'' and
			c.CONCEPT_PATH in ' + @parsed_exposure_concept_paths + 
			case when @patientset_instance_id is not null then 
		    ' and f.PATIENT_NUM in (select patient_num from QT_PATIENT_SET_COLLECTION 
		       where RESULT_INSTANCE_ID = ' + cast(@patientset_instance_id as varchar(max)) + ')' else ' ' end+
	' group by patient_num, isnull(c.concept_cd, left(c.name_char,50)), cast(start_date as DATE)
	union
	select patient_num, ''O'' concept_cd, cast(start_date as DATE) concept_date, row_number() over (partition by patient_num order by cast(start_date as DATE)) concept_order
	from observation_fact f, concept_dimension c, concept_dimension c2
	where f.CONCEPT_CD = c2.CONCEPT_CD and c2.CONCEPT_PATH like c.CONCEPT_PATH + ''%'' and
			c.CONCEPT_PATH in ' + @parsed_outcome_concept_paths +   
			case when @patientset_instance_id is not null then 
		    ' and f.PATIENT_NUM in (select patient_num from QT_PATIENT_SET_COLLECTION 
		       where RESULT_INSTANCE_ID = ' + cast(@patientset_instance_id as varchar(max)) + ')' else ' ' end+
	'group by patient_num, cast(start_date as DATE)
		union
	select patient_num, ''V'' concept_cd, cast(start_date as DATE) concept_date, row_number() over (partition by patient_num order by cast(start_date as DATE)) concept_order
	from observation_fact f, concept_dimension c, concept_dimension c2
	where f.CONCEPT_CD = c2.CONCEPT_CD and c2.CONCEPT_PATH like c.CONCEPT_PATH + ''%'' and
			c.CONCEPT_PATH in ' + @parsed_visit_concept_paths + 
			case when @patientset_instance_id is not null then 
		    ' and f.PATIENT_NUM in (select patient_num from QT_PATIENT_SET_COLLECTION 
		       where RESULT_INSTANCE_ID = ' + cast(@patientset_instance_id as varchar(max)) + ')' else ' ' end+
	'group by patient_num, cast(start_date as DATE)
	order by patient_num, concept_cd, concept_date
	'

	EXEC @Rcd = sp_executesql @exec_string	  

	IF @Rcd <> 0 OR @Rcd is NULL
	begin
		print 'Failed in concept_path metadata counts'
		return isnull (@Rcd, 99) -- return error code to the calling program
	end	
GO


